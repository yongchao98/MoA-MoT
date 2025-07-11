import math

# Part 1: Number of distinct planar graphs

# Plan for Part 1:
# 1. The problem asks for the number of distinct planar (non-crossing) graphs for 4-point scattering at 3 loops in phi^3 theory, excluding those with vertex corrections.
# 2. This standardly refers to "primitive" or "skeleton" graphs.
# 3. The number of such graphs at L-loop order is given by the L-th Catalan number, C_L.
# 4. The formula for the Catalan number is C_L = (1 / (L+1)) * binomial(2L, L).
# 5. We are at 3-loop order, so we will calculate C_3.

print("--- Part 1: Number of Planar Graphs ---")

L1 = 3  # Number of loops for Part 1

# Calculate binomial coefficient (2L choose L)
comb = math.comb(2 * L1, L1)

# Calculate the Catalan number C_L
num_diagrams = comb // (L1 + 1)

print(f"The number of diagrams is given by the Catalan number C_L for L={L1}.")
print(f"The formula is C_L = (1 / (L+1)) * (2L choose L).")
print(f"For L={L1}, the calculation is C_{L1} = (1 / ({L1}+1)) * ({2*L1} choose {L1})")
print(f"The binomial coefficient (6 choose 3) is {comb}.")
print(f"So, C_3 = (1/4) * {comb} = {num_diagrams}")
print(f"Result: The number of distinct planar graphs is {num_diagrams}.\n")


# Part 2: Power of the leading divergent term

# Plan for Part 2:
# 1. The problem asks for the power of the leading divergent term in the epsilon expansion of the Feynman integral.
# 2. For massless on-shell amplitudes, the dominant divergences are of Infrared (IR) origin.
# 3. The leading IR pole for an L-loop massless amplitude is of the order 1/epsilon^(2*L).
# 4. We are at 3-loop order, so L=3. The leading term is proportional to 1/epsilon^(2*3) = 1/epsilon^6.
# 5. A term like 1/epsilon^6 is written as epsilon^(-6), so the power is -6.

print("--- Part 2: Power of Leading Divergence ---")

L2 = 3  # Number of loops for Part 2

# The order of the leading IR pole in epsilon is 2*L
pole_order = 2 * L2

# The power of the leading divergent term is the exponent of epsilon, which is negative of the pole order.
power = -pole_order

print(f"For an L-loop massless on-shell amplitude, the leading IR divergence has a pole of order 2*L.")
print(f"For L = {L2} loops, the pole order is 2 * {L2} = {pole_order}.")
print(f"The leading term is proportional to 1/epsilon^{pole_order}, which is epsilon^({-pole_order}).")
print(f"Result: The power of the leading divergent term is {power}.")
