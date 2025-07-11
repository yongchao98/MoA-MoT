import math

# The problem reduces to finding r0 > 15 that satisfies the condition
# (T2(r0))^4 = 1/4, where T2(r) = (3*r - 37) / (r + 4).
# This gives |T2(r0)| = (1/4)**(1/4) = 1/sqrt(2).
# For r0 > 15, T2(r0) is positive, so T2(r0) = 1/sqrt(2).
# (3*r0 - 37) / (r0 + 4) = 1/sqrt(2)
# Let c = 1/sqrt(2).
# 3*r0 - 37 = c*(r0 + 4) = c*r0 + 4*c
# (3 - c)*r0 = 37 + 4*c
# r0 = (37 + 4*c) / (3 - c)

c = (1/4)**(1/4)
r0 = (37 + 4 * c) / (3 - c)

print(f"The radial distance r0 is calculated from the equation:")
print(f"(3 * r0 - 37) / (r0 + 4) = (1/4)**(1/4)")
print("\nLet c = (1/4)**(1/4). The equation for r0 is:")
print(f"(3 - c) * r0 = 37 + 4 * c")
print("\nThe numbers in this equation are:")
print(f"3, c={c}, 37, 4")
print(f"\nSolving for r0:")
print(f"r0 = (37 + 4 * {c}) / (3 - {c})")
print(f"r0 = {r0}")

# Final Answer
final_answer = r0