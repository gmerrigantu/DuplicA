# Use an official R base image
FROM rocker/r-ver:4.4.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    nodejs \
    npm \
    diamond-aligner \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    pkg-config \
    libfontconfig1-dev \
    libsodium-dev \
    automake \
    libxml2-dev \
    build-essential \
    git



# Install devtools and github packages
RUN R -e "install.packages('devtools', repos='https://cloud.r-project.org'); \
           devtools::install_gitlab('sandve-lab/evemodel'); \
           devtools::install_github('peterbchi/CNVSelectR')"


# Install CRAN packages
RUN R -e "install.packages('plumber', repos = 'https://cloud.r-project.org')"
RUN Rscript -e "stopifnot(requireNamespace('plumber', quietly = TRUE))"

RUN R -e "install.packages(c( \
  'bslib', 'htmltools', 'fs', 'tools', \
  'testthat', 'tidyverse', 'readxl', 'ape', 'R.utils', \
  'seqinr', 'igraph', 'biomartr', 'data.table', 'RColorBrewer', \
  'jsonlite', 'future', 'furrr' \
), repos = 'https://cloud.r-project.org')"


# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c('Biostrings', 'biomaRt', 'biomartr', 'clusterProfiler'), ask = FALSE)"


WORKDIR /app
COPY app /app

# Prefer system-provided diamond to avoid CPU instruction issues on emulation
RUN mv /app/Dependencies/OrthoFinder/bin/diamond /app/Dependencies/OrthoFinder/bin/diamond_avx || true \
 && ln -sf /usr/bin/diamond /app/Dependencies/OrthoFinder/bin/diamond

# Make OrthoFinder binaries executable
RUN chmod +x /app/Dependencies/OrthoFinder/orthofinder /app/Dependencies/OrthoFinder/bin/*

WORKDIR /app/frontend/duplic-a
RUN npm install && npm run build

RUN ls -la /app/frontend/duplic-a/public 




# Go back to app root
WORKDIR /app

# Expose ports
EXPOSE 8000 8001 8002

# Start all services
CMD ["sh", "-c", "\
  Rscript run_api.R & \
  python3 status_file_hosting.py & \
  cd frontend/duplic-a && npx serve -s public -l 8000 --no-clipboard & \
  wait"]

