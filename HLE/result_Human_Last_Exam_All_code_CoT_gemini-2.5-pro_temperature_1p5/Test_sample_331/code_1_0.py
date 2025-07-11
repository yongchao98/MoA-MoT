# Based on the properties of the Chow group of 0-cycles for each variety,
# the values for m(X) and M(X) are determined as follows:

# X1: Genus 2 curve
# The existence degree depends on whether the point is a Weierstrass point.
# m(X1) corresponds to edeg at a Weierstrass point.
# M(X1) corresponds to edeg at a general point.
m_X1 = 2
M_X1 = 3

# X2: General genus 7 curve
# For a general curve of genus g != 3, 5, the existence degree is g for all points.
m_X2 = 7
M_X2 = 7

# X3: Enriques surface
# An Enriques surface has a trivial group of 0-cycles of degree zero,
# which implies the existence degree is 1 for all points.
m_X3 = 1
M_X3 = 1

# X4: Grassmannian G(3,6)
# A Grassmannian is a rationally connected variety, which also implies
# its group of 0-cycles of degree zero is trivial.
# Thus, the existence degree is 1 for all points.
m_X4 = 1
M_X4 = 1

# Printing the final list of pairs as requested in the format.
# Each number is outputted by the print function via the formatted string.
print(f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})")
