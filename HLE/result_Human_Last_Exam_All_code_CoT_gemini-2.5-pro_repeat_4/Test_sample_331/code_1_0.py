# This script calculates the (m(X), M(X)) pairs for four given varieties
# based on known results in algebraic geometry.

# Case 1: X_1, a genus 2 curve
# A genus 2 curve C is hyperelliptic and has g=2. It has 2g+2 = 6 Weierstrass points.
# For a Weierstrass point p, the image of the curve in its Jacobian is symmetric,
# which leads to edeg(C, p) = 2.
# For any other point, edeg(C, p) = 3, which is the maximum possible value g+1.
# Therefore, the minimum value m(X_1) is 2 and the supremum M(X_1) is 3.
m_X1 = 2
M_X1 = 3

# Case 2: X_2, a general genus 7 curve
# For a general (non-hyperelliptic, etc.) curve C of genus g, a result by C. Voisin
# states that edeg(C, p) = g for all points p. Here g = 7.
# Thus, m(X_2) and M(X_2) are both equal to 7.
m_X2 = 7
M_X2 = 7

# Case 3: X_3, an Enriques surface
# For an Enriques surface S, the group of 0-cycles of degree 0 modulo rational
# equivalence, A_0(S), is isomorphic to Z/2Z.
# A direct analysis of the definition of edeg(S, p) for a group of order 2
# shows that edeg(S, p) = 2 for any point p on S.
# Hence, m(X_3) = 2 and M(X_3) = 2.
m_X3 = 2
M_X3 = 2

# Case 4: X_4, the Grassmannian G(3,6)
# The Grassmannian G(3,6) is a rational variety. For any smooth rational
# projective variety X, the Chow group of 0-cycles CH_0(X) is isomorphic to Z.
# This means any two points p, q are rationally equivalent, i.e., [p] = [q].
# From the definition, this implies edeg(X, p) = 1 for all points p.
# So, m(X_4) = 1 and M(X_4) = 1.
m_X4 = 1
M_X4 = 1

# Consolidate results into a list of pairs
results = [
    (m_X1, M_X1),
    (m_X2, M_X2),
    (m_X3, M_X3),
    (m_X4, M_X4)
]

# Format and print the final answer as specified.
# The f-string constructs the output from the computed numbers.
output_string = f"({results[0][0]}, {results[0][1]}), ({results[1][0]}, {results[1][1]}), ({results[2][0]}, {results[2][1]}), ({results[3][0]}, {results[3][1]})"
print(output_string)