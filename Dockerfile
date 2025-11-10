# Use Ubuntu as a base image (includes bash)
FROM ubuntu:22.04

# Prevent interactive prompts during package installs
ENV DEBIAN_FRONTEND=noninteractive

# Install bash and other dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash wget ca-certificates tar gzip && \
    rm -rf /var/lib/apt/lists/*

# Install SRA Toolkit
RUN mkdir -p /opt/sratoolkit && \
	wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz -C /opt/sratoolkit --strip-components=1 && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm sratoolkit.current-ubuntu64.tar.gz

# Default shell
CMD ["/bin/bash"]