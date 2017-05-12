Tutorial
========

1. Discuss HMMs - profile HMMs.
2. Discuss Pair-HMMs for alignment (S,M,D,I,E)
   T[s1,s2] is probability of transition from s1 -> s2 in R script.
3. Show paths:
   S(0,0) - M(1,1) - D(2,1) - D(3,1) - M(4,2) - I(4,3) - M(5,3) - E(5,3)
4. [note to self: determine notation for data, path]
   Maybe paths are pi[k] = (Z[k], I[k], J[k]) ?
   Maybe sequence data is X[1..L(X)] and Y[1..L(Y)]

5. Derive forward probability
   FM[i, j, z] = Pr(X_1:i[k], Y_1:j[k], Z[k] = z) for some k.
6. Show forward matrix, discuss slices M, D, I.

7. Derive backward probability
   Pr(S[k-1] | S[k], X, Y)
8. Sample a path and show an alignment.

9. Discuss why we want to sum over alignments when estimating delta, t

10. Discuss simple MCMC chain.
11. Run mcmc chain.