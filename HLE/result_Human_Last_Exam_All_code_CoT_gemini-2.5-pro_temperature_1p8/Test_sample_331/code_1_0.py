# Based on the mathematical reasoning outlined above, this script will now print the computed pairs.

# For X_1, a genus 2 curve: m(X_1) = 2, M(X_1) = 3
m_X1 = 2
M_X1 = 3

# For X_2, a general genus 7 curve: m(X_2) = 8, M(X_2) = 8
m_X2 = 8
M_X2 = 8

# For X_3, an Enriques surface: m(X_3) = 1, M(X_3) = 1
m_X3 = 1
M_X3 = 1

# For X_4, a Grassmannian G(3,6): m(X_4) = 1, M(X_4) = 1
m_X4 = 1
M_X4 = 1

# Format the final output string as requested.
# The format is (m(X_1), M(X_1)), (m(X_2), M(X_2)), (m(X_3), M(X_3)), (m(X_4), M(X_4))
final_answer = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"

print(final_answer)
