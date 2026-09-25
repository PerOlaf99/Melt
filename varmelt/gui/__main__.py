"""Entry point: ``python -m varmelt.gui [project.json]``."""
import sys

from .main import main

if __name__ == "__main__":
    sys.exit(main(sys.argv))