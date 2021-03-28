FROM rootproject/root:6.22.08-ubuntu20.04

WORKDIR /specfit/build

COPY ./ ../

ENV SPECFIT /specfit/build
ENV PATH ${SPECFIT}:${PATH}
ENV PYTHONPATH ${SPECFIT}:${PYTHONPATH}

RUN apt-get update -y
RUN apt-get install python3-pandas -y

RUN cmake ../\
    && make -j3

CMD ["./specfit.py", "--help"]
