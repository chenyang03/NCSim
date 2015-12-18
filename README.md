# NCSim

Simulator for Decentralized Network Coordinate Algorithms (NCSim) 
v1.0.0 

Author: Yang Chen (ychen@cs.duke.edu)


== Overview ==

Network Coordinate (NC) systems play an important role in scalable Internet distance prediction. NCSim is a MATLAB software package for simulating different decentralized Network Coordinate (NC) algorithms. With this tool, you can easily evaluation the performance of selected NC systems using traces collected from the Internet, and build your own applications on top of it.

== Coverage ==

Systems: Vivaldi (random), Vivaldi (height), Vivaldi (TIV aware), Phoenix, DMF, IDES
Evaluation Metrics: Relative Error, Rank Accuracy
Datasets: PL (169 nodes from Planetlab testbed), Toread (335 nodes from Planetlab testbed), King (1740 nodes)

Please execute 'NCSim_main.m' using MATLAB

== References ==

1. F. Dabek, R. Cox, and F. Kaashoek. Vivaldi: A Decentralized Network Coordinate System. In Proc. of ACM SIGCOMM, 2004.
2. G. Wang, B. Zhang, T.S.E. Ng. Towards Network Triangle Inequality Violation Aware Distributed Systems. In Proc. of ACM IMC, 2007.
3. Y. Chen, X. Wang, C. Shi, E.K. Lua, X. Fu, B. Deng, X. Li. Phoenix: A Weight-based Network Coordinate System Using Matrix Factorization. To appear in IEEE Transactions on Network and Service Management, 2011, 8(4).
4. Y. Liao, P. Geurts and G. Leduc. Network Distance Prediction Based on Decentralized Matrix Factorization. Proc. of IFIP Networking 2010.
5. Y. Mao, L. Saul, and J. M. Smith. IDES: An Internet Distance Estimation Service for Large Network. IEEE Journal on Selected Areas in Communications (JSAC), 2006.

== Note ==

It is difficult to conduct a fair comparison among different NC systems using the same framework, because different systems adopt different protocols for peer discovery, different considerations for handling node churn and distance variation. Moreover, many of them have several branches. Focusing more on NC calculate algorithms, this tool is not aiming at an exact prototype implementation for listed NC systems (please check with corresponding authors of these systems for an official release).
