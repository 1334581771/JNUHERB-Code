#!/usr/bin/env python
"""JNUHERB — Web application entry point.

Can be invoked directly as a standalone script (e.g. by Streamlit
Cloud) or via ``python -m jnuherb``.

Usage::

    streamlit run app_entrypoint.py
    # or
    python app_entrypoint.py
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from importlib.resources import files
from pathlib import Path


def _get_app_script() -> Path:
    """Locate the bundled Streamlit app script."""
    return files("jnuherb.data").joinpath("streamlit_app.py")


def main() -> None:
    """CLI entry point — parse args and launch the Streamlit server."""
    parser = argparse.ArgumentParser(
        description="JNUHERB — TCM Network Pharmacology Web App",
    )
    parser.add_argument(
        "--port", type=int, default=8501,
        help="TCP port number (default: 8501)",
    )
    parser.add_argument(
        "--host", type=str, default="127.0.0.1",
        help="Bind address. Use 0.0.0.0 for network access.",
    )
    parser.add_argument(
        "--no-browser", action="store_true",
        help="Do not automatically open a browser tab.",
    )
    parser.add_argument(
        "--theme-base", type=str, default="light",
        choices=["light", "dark"],
        help="Streamlit theme base.",
    )
    args = parser.parse_args()

    app_script = _get_app_script()
    if not app_script.is_file():
        raise FileNotFoundError(
            "Streamlit app not found in the jnuherb package. "
            "Reinstall the package or check that streamlit_app.py "
            "is bundled inside jnuherb/data/."
        )

    cmd = [
        sys.executable, "-m", "streamlit", "run",
        str(app_script),
        "--server.port", str(args.port),
        "--server.address", args.host,
        "--theme.base", args.theme_base,
    ]
    if args.no_browser:
        cmd.append("--server.headless=true")

    print(f"🌿 Starting JNUHERB Streamlit app...")
    print(f"   URL: http://{args.host}:{args.port}")
    print(f"   Press Ctrl+C to stop.")
    subprocess.run(cmd)


if __name__ == "__main__":
    main()
