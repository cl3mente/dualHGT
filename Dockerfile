# syntax=docker/dockerfile:1

# Comments are provided throughout this file to help you get started.
# If you need more help, visit the Dockerfile reference guide at
# https://docs.docker.com/go/dockerfile-reference/

# Want to help us make this template better? Share your feedback here: https://forms.gle/ybq9Krt8jtBL3iCk7

ARG PYTHON_VERSION=3.10.2
FROM python:${PYTHON_VERSION}-slim as base

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

# Download dependencies as a separate step to take advantage of Docker's caching.
# Leverage a cache mount to /root/.cache/pip to speed up subsequent builds.
# Leverage a bind mount to requirements.txt to avoid having to copy them into
# into this layer.
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=requirements.txt,target=requirements.txt \
    python -m pip install -r requirements.txt

# Include OrthoFinder in the image
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    mcl \
    fastme \
    diamond \
    wget \
    unzip \
    wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz && \
    tar xf OrthoFinder.tar.gz && \
    rm OrthoFinder.tar.gz

# Include KaKs_Calculator in the image

RUN wget https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download && \
    tar -xf KaKs_Calculator2.0.tar.gz && \
    rm KaKs_Calculator2.0.tar.gz

# Include paraAT in the image
RUN wget https://download.cncb.ac.cn/bigd/tools/ParaAT2.0.tar.gz && \
    tar -xf ParaAT2.0.tar.gz && \
    rm ParaAT2.0.tar.gz

# Switch to the non-privileged user to run the application.
USER appuser

# Copy the source code into the container.
COPY . .

# Expose the port that the application listens on.
EXPOSE 8000

# Run the application.
# CMD python HGT.py -i /data/bioinf2023/PlantPath2023/genomeANDgff -OFr /data/bioinf2023/PlantPath2023/genomeANDgff/results/prot/f7b812/Results_Feb23  -v -nt 50
