#
# This script calculates the pairs (m(X), M(X)) for four given algebraic varieties.
# The calculation is based on known results in algebraic geometry regarding the
# Chow group of 0-cycles for these varieties.
#

# Case 1: X_1, a genus 2 curve.
# For a genus 2 curve C, edeg(C, p) = 2 if p is one of the 6 Weierstrass points,
# and edeg(C, p) = 3 for any other point.
# Therefore, the minimum value m(X_1) is 2, and the maximum value M(X_1) is 3.
m_X1 = 2
M_X1 = 3

# Case 2: X_2, a general genus 7 curve.
# For a generic algebraic curve C of genus g >= 4, it is known that edeg(C, p) = g - 1
# for all points p on C. For g=7, this value is 7 - 1 = 6.
# Thus, m(X_2) and M(X_2) are both 6.
g_X2 = 7
m_X2 = g_X2 - 1
M_X2 = g_X2 - 1

# Case 3: X_3, an Enriques surface.
# For an Enriques surface S, the group of 0-cycles of degree 0 modulo rational
# equivalence is trivial (A_0(S) = 0). This implies that [p] = [q] in CH_0(S) for any
# two points p, q. So, for d=1, the equation 1*[p] = [q_1] holds for any p and q_1.
# Therefore, edeg(S, p) = 1 for all p, so m(X_3) = M(X_3) = 1.
m_X3 = 1
M_X3 = 1

# Case 4: X_4, the Grassmannian G(3,6).
# A Grassmannian is a rationally connected variety. For any rationally connected
# variety X, the Chow group of 0-cycles CH_0(X) is isomorphic to Z. This again
# implies [p] = [q] for any two points p, q.
# Similar to the Enriques surface, edeg(X, p) = 1 for all p.
# So, m(X_4) = M(X_4) = 1.
m_X4 = 1
M_X4 = 1

# Format and print the final list of pairs.
result_string = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"
print(result_string)