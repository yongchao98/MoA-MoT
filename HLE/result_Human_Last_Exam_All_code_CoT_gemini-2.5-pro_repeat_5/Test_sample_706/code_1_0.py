import math

# We are asked to find the limit of v(c) as c -> infinity.
# Based on the step-by-step derivation, the final answer is a constant.

# Let's outline the components of the final speed v(c):
# v(c) = pi_0(c) * v_0(c) + pi_1(c) * v_1(c)

# 1. Speed on the lower level v_0(c)
# v_0(c) = tanh(c). The limit as c -> inf is 1.
v0_inf = 1

# 2. Speed on the upper level v_1(c)
# In the limit c -> inf, this is the expected displacement per step.
# P(right edge exists) = 2/3. Displacement = +1.
# P(right edge missing) = 1/3.
#   P(vertical edge exists | right missing) = 1/2. Displacement = 0.
#   P(vertical edge missing | right missing) = 1/2. Displacement = -1.
v1_inf = (2/3) * 1 + (1/3) * (1/2) * 0 + (1/3) * (1/2) * (-1)
# v1_inf = 2/3 - 1/6 = 3/6 = 1/2

# 3. Stationary probabilities pi_0(c) and pi_1(c)
# The ratio pi_1(c)/pi_0(c) is proportional to exp(-c).
# So, as c -> inf, pi_1(c) -> 0 and pi_0(c) -> 1.
pi0_inf = 1
pi1_inf = 0

# 4. Final calculation
# lim v(c) = (lim pi_0(c)) * (lim v_0(c)) + (lim pi_1(c)) * (lim v_1(c))
final_speed = pi0_inf * v0_inf + pi1_inf * v1_inf

# The code will now print the final equation with the computed values.
# Note that the values are for the limit c -> infinity.
print("The asymptotic speed v is calculated as the limit of v(c) as c approaches infinity.")
print("v(c) = π₀(c) * v₀(c) + π₁(c) * v₁(c)")
print("\nAs c → ∞, the limits of the components are:")
print(f"lim π₀(c) = {pi0_inf}")
print(f"lim v₀(c) = {v0_inf} (since lim tanh(c) = 1)")
print(f"lim π₁(c) = {pi1_inf}")
print(f"lim v₁(c) = 2/3 - 1/6 = {v1_inf}")
print("\nSubstituting these limits into the equation:")
print(f"v = ({pi0_inf}) * ({v0_inf}) + ({pi1_inf}) * ({v1_inf})")
print(f"v = {pi0_inf * v0_inf} + {pi1_inf * v1_inf}")
print(f"v = {final_speed}")
