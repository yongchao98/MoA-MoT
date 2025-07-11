# My thinking process to arrive at the values is as follows:
# The value edeg(X,p) is the smallest integer d such that for any q in X,
# the cycle d[p] - [q] is equivalent to an effective cycle of degree d-1.

# Case X_1: A genus g=2 curve C.
# The value edeg(C,p) is the smallest d such that h^0(dp-q) >= 1 for all q in C.
# By Riemann-Roch, h^0(dp-q) >= deg(dp-q) - g + 1 = d-1-g+1 = d-g.
# If d = g+1 = 3, h^0(3p-q) >= 3-2 = 1. So edeg(C,p) <= 3 for all p.
# If d = g = 2, we need h^0(2p-q) >= 1 for all q.
# h^0(2p-q) = h^0(K_C-2p+q).
# If p is a Weierstrass point (6 of them on a g=2 curve), 2p is the canonical divisor K_C.
# Then h^0(K_C-2p+q) = h^0(q) = 1. So edeg(C,p) = 2 for Weierstrass points.
# If p is not a Weierstrass point, then for a general q, h^0(K_C-2p+q) = 0.
# So edeg(C,p) > 2, which means edeg(C,p) = 3.
# Therefore, the minimum value m(X_1) is 2 and the maximum M(X_1) is 3.
g1 = 2
m1 = g1
M1 = g1 + 1

# Case X_2: A general genus g=7 curve C.
# For a general curve of genus g, it is a known result that edeg(C,p) = g+1 for all points p.
# This is because for a general curve and general point p, the required positivity
# of h^0 fails for d=g.
g2 = 7
m2 = g2 + 1
M2 = g2 + 1

# Case X_3: An Enriques surface S.
# For an Enriques surface, the Chow group of 0-cycles of degree zero, A_0(S), is trivial.
# This means [p] is rationally equivalent to [q] for any points p,q.
# The condition for edeg(S,p)=d is that for any q_1, there exist q_2,...,q_d such that
# d[p] is rationally equivalent to q_1+...+q_d.
# Since [q_i] ~ [p] for all i, this becomes d[p] ~ d[p], which is always true.
# We need the minimum positive d, which is 1.
# For d=1, the condition is [p] ~ [q_1] for all q_1, which is true since A_0(S)=0.
# So edeg(S,p) = 1 for all p.
m3 = 1
M3 = 1

# Case X_4: Grassmannian G(3,6).
# A Grassmannian is a rationally connected variety. For any rationally connected
# variety X, the Chow group A_0(X) is trivial.
# The argument is identical to the Enriques surface case.
# edeg(X,p) = 1 for all p.
m4 = 1
M4 = 1

# Printing the final result in the specified format.
print(f"({m1}, {M1}), ({m2}, {M2}), ({m3}, {M3}), ({m4}, {M4})")