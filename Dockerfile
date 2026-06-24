# =============================================================================
# JNUHERB — Network Pharmacology Web Application (Docker)
# =============================================================================
# Build:
#   docker build -t jnuherb .
#
# Run:
#   docker run -p 8501:8501 jnuherb
#
# With persistent results:
#   docker run -p 8501:8501 -v $(pwd)/results:/home/jnuherb/results jnuherb
# =============================================================================

FROM python:3.11-slim

LABEL org.opencontainers.image.title="jnuherb"
LABEL org.opencontainers.image.description="TCM Network Pharmacology Analysis Web App"
LABEL org.opencontainers.image.version="0.1.0"

# ---- System dependencies ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# ---- App user (non-root) ----
RUN useradd --create-home --shell /bin/bash jnuherb
USER jnuherb
WORKDIR /home/jnuherb

ENV PATH="/home/jnuherb/.local/bin:${PATH}"

# ---- Install Python dependencies (layer caching) ----
COPY --chown=jnuherb:jnuherb pyproject.toml /home/jnuherb/app/pyproject.toml
RUN pip install --no-cache-dir --user \
    pandas>=1.5.0 \
    numpy>=1.24.0 \
    matplotlib>=3.7.0 \
    seaborn>=0.12.0 \
    plotly>=5.14.0 \
    streamlit>=1.28.0

# ---- Copy and install the package ----
COPY --chown=jnuherb:jnuherb src/ /home/jnuherb/app/src/
COPY --chown=jnuherb:jnuherb README.md LICENSE /home/jnuherb/app/
RUN pip install --no-cache-dir --user /home/jnuherb/app/

# ---- Streamlit config ----
RUN mkdir -p /home/jnuherb/.streamlit
COPY --chown=jnuherb:jnuherb .streamlit/config.toml /home/jnuherb/.streamlit/config.toml

# ---- Results volume ----
RUN mkdir -p /home/jnuherb/results
VOLUME /home/jnuherb/results

# ---- Health check ----
HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 \
    CMD curl -sf http://localhost:8501/_stcore/health || exit 1

EXPOSE 8501

# ---- Launch via CLI entry point ----
ENTRYPOINT ["jnuherb", "--host", "0.0.0.0", "--no-browser"]
CMD ["--port", "8501"]
