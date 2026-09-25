"""Shared HTTP helper: polite, retrying JSON fetches.

UCSC and NCBI throttle the API.  A global token-bucket rate limiter keeps
batches (which fetch in parallel threads) well under the limits, and every
request retries with exponential backoff on ``429 Too Many Requests`` and
transient 5xx/network errors so variants finish instead of failing.
"""

import json
import threading
import time
import urllib.error
import urllib.request

RATE_PER_SECOND = 3.0     # max sustained requests per second (UCSC/NCBI polite)
BURST = 6                 # short-term burst allowed
MAX_RETRIES = 6

USER_AGENT = "varmelt/1.0"


class _RateLimiter:
    """Thread-safe token bucket shared by all HTTP clients in this process."""

    def __init__(self, rate: float, burst: float):
        self._rate = rate
        self._burst = float(burst)
        self._tokens = float(burst)
        self._last = time.monotonic()
        self._lock = threading.Lock()

    def acquire(self):
        while True:
            with self._lock:
                now = time.monotonic()
                self._tokens = min(self._burst,
                                   self._tokens + (now - self._last) * self._rate)
                self._last = now
                if self._tokens >= 1.0:
                    self._tokens -= 1.0
                    return
            time.sleep(0.01)


_LIMITER = _RateLimiter(RATE_PER_SECOND, BURST)


def http_get_json(url: str, headers: dict = None, timeout: int = 20,
                  retries: int = MAX_RETRIES) -> dict:
    """Fetch *url* and parse its JSON body, retrying on throttling errors.

    Re-raises the last exception if every attempt fails.
    """
    hdr = {"User-Agent": USER_AGENT}
    if headers:
        hdr.update(headers)
    req = urllib.request.Request(url, headers=hdr)

    attempt = 0
    while True:
        _LIMITER.acquire()
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as exc:
            if exc.code != 429 and exc.code < 500:
                raise
            attempt += 1
            if attempt >= retries:
                raise
            time.sleep(_backoff(attempt, exc.headers.get("Retry-After")))
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            attempt += 1
            if attempt >= retries:
                raise
            time.sleep(_backoff(attempt, None))


def http_post(url: str, data: bytes, headers: dict = None, timeout: int = 30,
              retries: int = MAX_RETRIES) -> bytes:
    """POST *data* (form-encoded bytes) and return the response bytes.

    Retries on ``429``/5xx/network errors like the GET helper.
    """
    hdr = {"User-Agent": USER_AGENT, "Content-Type":
           "application/x-www-form-urlencoded"}
    if headers:
        hdr.update(headers)
    req = urllib.request.Request(url, data=data, headers=hdr)

    attempt = 0
    while True:
        _LIMITER.acquire()
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code != 429 and exc.code < 500:
                raise
            attempt += 1
            if attempt >= retries:
                raise
            time.sleep(_backoff(attempt, exc.headers.get("Retry-After")))
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            attempt += 1
            if attempt >= retries:
                raise
            time.sleep(_backoff(attempt, None))


def _backoff(attempt: int, retry_after) -> float:
    """Seconds to wait: honour Retry-After (capped) else exponential backoff."""
    if retry_after is not None:
        try:
            return min(float(retry_after), 30.0)
        except (TypeError, ValueError):
            pass
    return min(2.0 ** attempt, 30.0)


def set_rate(rate: float):
    """Reconfigure the global request rate (calls/sec)."""
    _LIMITER._rate = max(0.05, float(rate))