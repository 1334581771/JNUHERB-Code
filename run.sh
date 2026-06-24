#!/usr/bin/env bash
# =============================================================================
# JNUHERB — Quick launch script
# =============================================================================
# Usage:
#   ./run.sh                  # launch on http://127.0.0.1:8501
#   ./run.sh --port 8080      # custom port
#   ./run.sh --public         # listen on 0.0.0.0 (accessible from LAN)
#   ./run.sh --docker         # build & run via Docker
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PORT=8501
HOST="127.0.0.1"
MODE="local"

# ---- Parse args ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --port)
            PORT="$2"
            shift 2
            ;;
        --public)
            HOST="0.0.0.0"
            shift
            ;;
        --docker)
            MODE="docker"
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [--port N] [--public] [--docker]"
            echo ""
            echo "  --port N     TCP port (default 8501)"
            echo "  --public     Listen on 0.0.0.0 (LAN accessible)"
            echo "  --docker     Build and run via Docker"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ---- Launch ----
case "$MODE" in
    docker)
        echo "🐳 Building Docker image..."
        docker build -t jnuherb .
        echo "🚀 Starting JNUHERB container on port $PORT..."
        docker run --rm -p "$PORT:8501" jnuherb
        ;;
    local)
        # Ensure package is installed
        if ! python -c "import jnuherb" &>/dev/null; then
            echo "📦 Installing jnuherb in development mode..."
            pip install -e . --quiet
        fi
        echo "🌿 Launching jnuherb..."
        echo "   URL: http://${HOST}:${PORT}"
        echo "   Press Ctrl+C to stop."
        python app_entrypoint.py --host "$HOST" --port "$PORT"
        ;;
esac
