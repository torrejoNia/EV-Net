FROM rocker/rstudio:4.5.3

RUN apt-get update && apt-get install -y \
    curl \
    git \
    pkg-config \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libcairo2-dev \
    libxt-dev \
    libuv1-dev \
    libglpk-dev \
    libhdf5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error renv yaml

RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/latest'))" \
    >> /usr/local/lib/R/etc/Rprofile.site \
    && echo "options(HTTPUserAgent = sprintf('R/%s R (%s)', getRversion(), paste(getRversion(), R.version\$platform, R.version\$arch, R.version\$os)))" \
    >> /usr/local/lib/R/etc/Rprofile.site

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/
COPY --chown=rstudio:rstudio docker-renv /home/rstudio

RUN echo 'source("renv/activate.R")' > /home/rstudio/.Rprofile \
    && chown rstudio:rstudio /home/rstudio/.Rprofile

USER rstudio
RUN R -e "setwd('/home/rstudio'); source('renv/activate.R'); renv::restore(prompt = FALSE); renv::install('.')"
USER root