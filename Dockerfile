# Use bioconductor container as a base image
FROM bioconductor/bioconductor_docker:devel

WORKDIR /cromwell_root

# Prevent interactive prompts during package installs
ENV DEBIAN_FRONTEND=noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN=true

RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        g++ \
        zip \
        unzip \
        make \
        libxml2-dev \
        libpng-dev \
        libjpeg-dev \
        libcurl4-openssl-dev \
        zlib1g-dev \
        libbz2-dev \
        libpcre3-dev \
        libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Install desired Bioconductor packages
RUN R -e "BiocManager::install(c('ensembldb', 'EnsDb.Hsapiens.v75', 'tximport', 'edgeR', 'goseq', 'GO.db', 'org.Hs.eg.db', 'gplots', 'ggplot2', 'ggrepel', 'ghibli', 'GenomicFeatures'), ask = FALSE, update = TRUE)"

RUN rm -rf /tmp/Rtmp*
