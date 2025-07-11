# The user wants to solve a physics problem using the Titan computer architecture.
# The primary task is to determine if the calculation is feasible under the given constraints.

# 1. Physical Model Derivation
# Net force F_net = (F*cos(45), F*sin(45) - m*g)
# Acceleration a = F_net / m = (F*cos(45)/m, F*sin(45)/m - g)
# Trajectory: x = 0.5*a_x*t^2, y = 0.5*a_y*t^2
# Ratio: y/x = a_y / a_x = (F*sin(45)/m - g) / (F*cos(45)/m)
# Since sin(45)=cos(45), let's call it c.
# y/x = (F*c/m - g) / (F*c/m) = 1 - (m*g)/(F*c)
# (m*g)/(F*c) = 1 - y/x
# F = (m*g) / (c * (1 - y/x))
# With x=20, y=10, 1-y/x = 0.5. c = cos(45) = 1/sqrt(2)
# F = (m*g) / ((1/sqrt(2)) * 0.5) = 2*m*g*sqrt(2)

# Mass calculation:
# r = 0.5 cm = 0.005 m
# rho = 0.9 kg/cm^3 = 900,000 kg/m^3
# V = (4/3)*pi*r^3 = (4/3)*pi*(0.005)^3
# m = rho * V = 900000 * (4/3) * pi * (0.005)^3 = 0.15 * pi kg (~0.4712 kg)
# F = 2 * (0.15*pi) * g * sqrt(2) = 0.3 * pi * g * sqrt(2)

# 2. Feasibility Analysis under Titan Rules
# We need to calculate F = 0.3 * pi * g * sqrt(2) using fractions.
# Let's choose some reasonable approximations:
# pi ~ 22/7 (~3.14)
# g ~ 29/3 (~9.67)
# sqrt(2) ~ 24/17 (~1.41)
# 0.3 = 3/10
#
# The formula is F = (3/10) * (22/7) * (29/3) * (24/17)
# Let's perform the calculation step-by-step, checking constraints.
# Step 1: (3/10) * (22/7) = 66/70.
#   - Numerator 66 and denominator 70 both exceed the 5-bit limit of 31.
#   - We can simplify 66/70 to 33/35.
#   - Numerator 33 and denominator 35 still exceed the limit of 31.
#
# Let's try another path. Let's calculate the term pi * g first.
# pi * g ~ (22/7) * (29/3) = 638/21.
#   - Numerator 638 vastly exceeds 31. The fraction is irreducible. This intermediate
#     value cannot be stored.
#
# Let's try to calculate the value of the term `pi * g`.
# It is approximately 3.14 * 9.8 = 30.772.
# Let's try to multiply this by another term, 0.3.
# The product is 30.772 * 0.3 = 9.23.
# Even if we could calculate pi*g, we would then multiply it by other terms,
# and the intermediate values would only grow.
#
# The core problem is that to calculate the force, we must multiply several numbers together.
# For example, let's just calculate `2 * g * m`, where m ~ 0.47.
# The value is roughly `2 * 9.8 * 0.47 = 9.2`. This seems manageable.
# But `m` itself must be calculated from `(4/3) * pi * (rho * r^3)`.
# Let's analyze `rho * r^3`. `r_cm = 1/2`, `rho_kg_cm3 = 9/10`.
# `(r_cm)^3 = 1/8`. `rho * r^3 = (9/10) * (1/8) = 9/80`. Denominator `80 > 31`. This fails.
#
# It is impossible to represent the input physical parameters or to calculate intermediate
# products without violating the 5-bit integer constraint. Any valid sequence of
# operations under the Titan rules will fail.
# Therefore, Titan cannot be used to calculate the required force.

print("<<<N0>>>")