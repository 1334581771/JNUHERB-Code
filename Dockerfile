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

# ---- Stage 1: build binary wheels for the docking stack ----
# vina is a C++/SWIG extension linked against boost; rdkit ships manylinux
# wheels but pulling everything through one builder keeps the final image
# free of compilers. Runtime boost libraries ARE still required (see stage 2).
FROM python:3.11-slim AS dockbuilder

# Same Aliyun Debian mirror as the runtime stage — this stage downloads the
# most (full boost dev headers), so it must not hit deb.debian.org directly.
RUN sed -i 's|http://deb.debian.org/debian|https://mirrors.aliyun.com/debian|g' /etc/apt/sources.list.d/debian.sources \
    && sed -i 's|http://deb.debian.org/debian-security|https://mirrors.aliyun.com/debian-security|g' /etc/apt/sources.list.d/debian.sources \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       g++ swig libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip wheel --no-cache-dir --wheel-dir /wheels \
    -i https://mirrors.aliyun.com/pypi/simple/ \
    "rdkit>=2023.3" "meeko>=0.5" "gemmi>=0.6" "scipy>=1.10" "prody>=2.4" \
    "py3Dmol>=2.0" "vina>=1.2.3"

# ---- Stage 2: runtime image ----
FROM python:3.11-slim

LABEL org.opencontainers.image.title="jnuherb"
LABEL org.opencontainers.image.description="TCM Network Pharmacology Analysis Web App"
LABEL org.opencontainers.image.version="0.1.0"

# ---- System dependencies ----
RUN sed -i 's|http://deb.debian.org/debian|https://mirrors.aliyun.com/debian|g' /etc/apt/sources.list.d/debian.sources \
    && sed -i 's|http://deb.debian.org/debian-security|https://mirrors.aliyun.com/debian-security|g' /etc/apt/sources.list.d/debian.sources \
    && apt-get update \
    && apt-get install -y --no-install-recommends curl fontconfig \
       # runtime boost libraries needed by the vina extension module
       # (must match the libboost-all-dev default in the builder stage)
       libboost-thread1.88.0 libboost-serialization1.88.0 \
       libboost-filesystem1.88.0 libboost-program-options1.88.0 \
    && rm -rf /var/lib/apt/lists/*

# ---- Install bundled Chinese font so PDF exports can render CJK characters ----
COPY src/jnuherb/data/fonts/NotoSansCJKsc-Regular.otf /usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf
RUN fc-cache -f -v

# ---- App user (non-root) ----
RUN useradd --create-home --shell /bin/bash jnuherb
USER jnuherb
WORKDIR /home/jnuherb

ENV PATH="/home/jnuherb/.local/bin:${PATH}"

# Configure pip to use Aliyun mirror for faster builds inside China
ENV PIP_INDEX_URL=https://mirrors.aliyun.com/pypi/simple/
ENV PIP_TRUSTED_HOST=mirrors.aliyun.com
ENV PIP_TIMEOUT=120
ENV PIP_DEFAULT_TIMEOUT=120

# ---- Install Python dependencies (layer caching) ----
COPY --chown=jnuherb:jnuherb pyproject.toml /home/jnuherb/app/pyproject.toml
RUN pip install --no-cache-dir --user \
    pandas>=1.5.0 \
    numpy>=1.24.0 \
    matplotlib>=3.7.0 \
    seaborn>=0.12.0 \
    plotly>=5.14.0 \
    # pinned to the dev-tested version; the app relies on st.navigation,
    # st.logo and st.page_link (>=1.36/1.39) and AppTest-verified behaviour
    streamlit==1.58.0 \
    gseapy>=1.0.0 \
    requests>=2.28.0 \
    kaleido==0.2.1

# ---- Docking stack (pre-built wheels from stage 1) ----
COPY --from=dockbuilder --chown=jnuherb:jnuherb /wheels /tmp/wheels
RUN pip install --no-cache-dir --user /tmp/wheels/*.whl && rm -rf /tmp/wheels

# ---- Copy and install the package ----
COPY --chown=jnuherb:jnuherb src/ /home/jnuherb/app/src/
COPY --chown=jnuherb:jnuherb README.md LICENSE /home/jnuherb/app/
RUN pip install --no-cache-dir --user /home/jnuherb/app/

# ---- Bundled SQLite database ----
# Build jnuherb.db from the CSV tables inside the installed package, then
# drop the CSV copies (the app reads from the db; ensembl_uniprot.csv stays
# because the docking module loads it directly).
RUN python -m jnuherb.build_db && \
    find /home/jnuherb/.local -path '*/jnuherb/data/jnuherb_*.csv' ! -name 'ensembl_uniprot.csv' -delete && \
    rm -f /home/jnuherb/app/src/jnuherb/data/jnuherb_*.csv && \
    ls -lh /home/jnuherb/.local/lib/python3.11/site-packages/jnuherb/data/jnuherb.db

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
