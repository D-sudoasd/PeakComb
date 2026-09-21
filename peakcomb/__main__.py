"""python -m peakcomb launches the GUI; python -m peakcomb.cli runs headless."""

from __future__ import annotations


def main() -> None:
    from .gui.app import main as gui_main

    gui_main()


if __name__ == "__main__":
    main()
