# NCSim

Simulator for Decentralized Network Coordinate Algorithms (NCSim) 

Version 1.1.0 
Updated on Jan. 3, 2016

Copyright (C) <2011-2016> by Yang Chen, Fudan University (chenyang@fudan.edu.cn)


## Overview

Network Coordinate (NC) systems play an important role in scalable Internet distance prediction. NCSim is a MATLAB software package for simulating different decentralized Network Coordinate (NC) algorithms. With this tool, you can easily evaluation the performance of selected NC systems using traces collected from the Internet, and build your own applications on top of it.

## Coverage

Systems: Vivaldi (random), Vivaldi (height), Vivaldi (TIV aware), Phoenix, DMFSGD, IDES
Evaluation Metrics: Relative Error, Rank Accuracy
Datasets: PL (169 nodes from Planetlab testbed), Toread (335 nodes from Planetlab testbed), King (1740 nodes)

Please run 'src/NCSim_main.m' using MATLAB

## References

1. F. Dabek, R. Cox, and F. Kaashoek. Vivaldi: A Decentralized Network Coordinate System. In Proc. of ACM SIGCOMM, 2004.
2. G. Wang, B. Zhang, T.S.E. Ng. Towards Network Triangle Inequality Violation Aware Distributed Systems. In Proc. of ACM IMC, 2007.
3. Y. Chen, X. Wang, C. Shi, E.K. Lua, X. Fu, B. Deng, X. Li. Phoenix: A Weight-based Network Coordinate System Using Matrix Factorization. IEEE Transactions on Network and Service Management, 2011, 8(4):334-347.
4. Yongjun Liao, Wei Du, Pierre Geurts and Guy Leduc. DMFSGD: A Decentralized Matrix Factorization Algorithm for Network Distance Prediction. IEEE/ACM Transactions on Networking, 2013, 21(5):1511-1524.
5. Y. Mao, L. Saul, and J. M. Smith. IDES: An Internet Distance Estimation Service for Large Network. IEEE Journal on Selected Areas in Communications (JSAC), 2006, 24(12):2273-2284.

## Note

This package includes my implementation for several representative NC systems. My goal is to help the users study these systems using MATLAB conveniently. If you are interested in the offical release of Vivaldi, IDES and DMFSGD, please contact their authors.


## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT). 
