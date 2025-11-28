FROM mambaorg/micromamba:latest

# Copy environment specification
COPY environment.yml /tmp/environment.yml

# Install your environment into base
RUN micromamba install -y -n base -f /tmp/environment.yml --channel conda-forge \
    && micromamba clean --all --yes

# Optional: add your code
WORKDIR /workspace
COPY . .

SHELL ["/usr/bin/bash", "-c"]

CMD ["bash"]
