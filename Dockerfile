FROM mambaorg/micromamba:latest

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Create the environment in the base environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set working directory
WORKDIR /work

# Use bash as the default shell
SHELL ["/usr/bin/bash", "-c"]
