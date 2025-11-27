FROM mambaorg/micromamba:latest

COPY environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /work
SHELL ["/usr/bin/bash", "-c"]
