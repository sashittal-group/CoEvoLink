FROM mambaorg/micromamba:latest

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Install environment
RUN micromamba install -y -n base \
      -f /tmp/environment.yml \
      --channel conda-forge \
      && micromamba clean --all --yes

WORKDIR /work
SHELL ["/usr/bin/bash", "-c"]
