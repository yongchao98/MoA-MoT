import math

# Define the parameters given in the problem
# Let m=1 and c=1 for simplicity in calculation. The final answer will be in terms of m and c.
m = 1
c = 1

# Astronaut 1
m1 = 2 * m
v1 = 0.6 * c

# Astronaut 2
m2 = 3 * m
v2 = 0.5 * c

# Step 1: Calculate the relative speed (v_rel)
# The formula for relativistic velocity subtraction is: v_rel = (v_a - v_b) / (1 - (v_a * v_b) / c^2)
# Since both are moving in the same direction, we can use their speeds directly.
v_rel = (v1 - v2) / (1 - (v1 * v2) / c**2)

# Step 2: Calculate the total energy (E_total)
# The formula for relativistic energy of a single particle is E = gamma * m_rest * c^2,
# where gamma = 1 / sqrt(1 - v^2/c^2).

# Calculate gamma and energy for Astronaut 1
gamma1 = 1 / math.sqrt(1 - v1**2 / c**2)
E1 = gamma1 * m1 * c**2

# Calculate gamma and energy for Astronaut 2
gamma2 = 1 / math.sqrt(1 - v2**2 / c**2)
E2 = gamma2 * m2 * c**2

# The total energy of the system is the sum of the individual energies.
E_total = E1 + E2

print(f"Relative speed: {v_rel:.2f}c")
print(f"Total energy: {E_total:.2f} mc^2")