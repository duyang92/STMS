# STMS

Short-term memory sampling (STMS) is a novel spread estimator. It will be presented at the IEEE International Conference on Computer Communications ([INFOCOM 2022](https://infocom2022.ieee-infocom.org/)). This is the repository for STMS's source codes.

## Introduction
Per-flow spread measurement in high-speed networks can provide indispensable information to many practical applications. However, it is challenging to measure millions of flows at line speed because on-chip memory modules cannot simultaneously provide large capacity and large bandwidth. The prior studies address this mismatch by entirely using on-chip compact data structures or utilizing off-chip space to assist limited on-chip memory. 
Nevertheless, their on-chip data structures record massive transient elements, each of which only appears in a short time interval in a long-period measurement task, and thus waste significant on-chip space. This paper presents short-term memory sampling, a novel spread estimator that samples new elements while only holding elements for short periods. Our estimator can work with tiny on-chip space and provide accurate estimations for online queries. The key of our design is a short-term memory duplicate filter that reports new elements and filters duplicates effectively while allowing incoming elements to override the stale elements to reduce on-chip memory usage. We implement our approach on a NetFPGA-equipped prototype. Experimental results based on real Internet traces show that, compared to the state-of-the-art, short-term memory sampling reduces up to 99% of on-chip memory usage when providing the same probabilistic assurance on spread-estimation error. 

## About this repo
### CPU
STMS is implemented on CPU, the code is compiled using gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04). To achieve a higher throughput, please turn on the O2 optimization option while compiling. The compiling command is as follows.
```shell
g++ STMS.cpp MurmurHash3.cpp -O2 -m64
```

### FPGA
STMS is implemented on FPGA.