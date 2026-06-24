"""Streamlit web application launcher for JNUHERB.

Can be used programmatically::

    >>> import jnuherb
    >>> jnuherb.launch_app(port=8501)

Or from the command line::

    python -m jnuherb --port 8501 --host 0.0.0.0
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from importlib.resources import files


def launch_app(
    port: int = 8501,
    host: str = "127.0.0.1",
    browser: bool = True,
) -> None:
    """Launch the JNUHERB Streamlit web application.

    Args:
        port: TCP port number (default 8501).
        host: Bind address.  Use ``"0.0.0.0"`` to allow network access.
        browser: Whether to automatically open a browser tab.
    """
    app_script = files("jnuherb.data").joinpath("streamlit_app.py")

    if not app_script.is_file():
        raise FileNotFoundError(
            "Streamlit app not found in the jnuherb package. "
            "The app may not be included in this version."
        )

    cmd = [
        sys.executable, "-m", "streamlit", "run",
        str(app_script),
        "--server.port", str(port),
        "--server.address", host,
    ]
    if not browser:
        cmd.append("--server.headless=true")

    print(f"Starting JNUHERB Streamlit app...")
    print(f"App URL: http://{host}:{port}")

    subprocess.run(cmd)


def main() -> None:
    """CLI entry point with argument parsing."""
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
    args = parser.parse_args()

    launch_app(port=args.port, host=args.host, browser=not args.no_browser)
