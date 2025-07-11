# Plan:
# 1. For X_1, a genus 2 curve: The curve is hyperelliptic.
#    - For the 6 Weierstrass points p, edeg(X_1, p) = 2.
#    - For any other point p, edeg(X_1, p) = 3.
#    - Thus, m(X_1) = 2 and M(X_1) = 3.
m_X1 = 2
M_X1 = 3

# 2. For X_2, a general genus 7 curve C:
#    - M(C) is the genus g, which is 7.
#    - m(C) is gon(C) + 1. For a general genus g=7 curve, gon(C) = floor((7+3)/2) = 5.
#    - So, m(X_2) = 5 + 1 = 6.
#    - Thus, m(X_2) = 6 and M(X_2) = 7.
m_X2 = 6
M_X2 = 7

# 3. For X_3, an Enriques surface S:
#    - The Chow group of 0-cycles of degree 0, A_0(S), is Z/2Z.
#    - This implies 2[p] is rationally equivalent to 2[q] for any points p, q.
#    - edeg(S, p) cannot be 1 as A_0(S) is not trivial.
#    - For d=2, the condition 2[p] = [q_1] + [q_2] is satisfied for any q_1 by choosing q_2 = q_1, since 2[p] ~ 2[q_1].
#    - This holds for all points p, so edeg(S, p) = 2 for all p.
#    - Thus, m(X_3) = 2 and M(X_3) = 2.
m_X3 = 2
M_X3 = 2

# 4. For X_4, the Grassmannian G(3,6):
#    - This is a rationally connected variety, so any two points p, q are rationally equivalent: [p] = [q].
#    - The condition for d=1 is: [p] = [q_1] for any q_1. This is true.
#    - So, edeg(X_4, p) = 1 for all p.
#    - Thus, m(X_4) = 1 and M(X_4) = 1.
m_X4 = 1
M_X4 = 1

# Format the final output string as requested.
# The format is (m(X_1), M(X_1)), (m(X_2), M(X_2)), (m(X_3), M(X_3)), (m(X_4), M(X_4))
final_answer = f"({m_X1}, {M_X1}), ({m_X2}, {M_X2}), ({m_X3}, {M_X3}), ({m_X4}, {M_X4})"

print(final_answer)