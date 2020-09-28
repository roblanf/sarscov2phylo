FROM continuumio/miniconda3:4.8.2

RUN conda update -n base -c defaults conda
COPY environment.yml /tmp/environment.yml

ARG USER_ID

RUN useradd -m --uid $USER_ID --shell /bin/bash user
RUN mkdir -p /work/sarscov2phylo
RUN chown user /work/sarscov2phylo
USER user

RUN conda env create -n sarscov2phylo -f /tmp/environment.yml

ADD . /work/sarscov2phylo
WORKDIR /work/sarscov2phylo
