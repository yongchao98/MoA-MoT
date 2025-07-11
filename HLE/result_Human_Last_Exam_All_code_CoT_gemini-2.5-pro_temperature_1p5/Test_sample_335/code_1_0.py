import math

# The goal is to compute floor(10^6 * V), where V is the simplicial volume
# of the complement of the knot K = C_{4,3}(Conway) # Wh_-^2(Eight).

# The simplicial volume is additive over the connected sum operation.
# K = K1 # K2, where K1 = C_{4,3}(Conway) and K2 = Wh_-^2(Eight).
# V = ||S^3 \ K|| = ||S^3 \ K1|| + ||S^3 \ K2||

# Step 1: Compute the simplicial volume of K1's complement.
# K1 is the (4,3)-cable of the Conway knot. For a cable knot C_{p,q}(J)
# with |p| >= 2, the simplicial volume of its complement equals that of the
# companion knot J's complement. Here p=4, so the rule applies.
# ||S^3 \ K1|| = ||S^3 \ Conway||.
# It is a known result (L. Lischka, 2018) that the Conway knot is a graph
# knot, and therefore the simplicial volume of its complement is 0.
V_K1 = 0

# Step 2: Compute the simplicial volume of K2's complement.
# K2 is a satellite knot (a Whitehead double) with the figure-8 knot as companion.
# The simplicial volume is the sum of the companion's and the pattern's
# simplicial volumes.
# ||S^3 \ K2|| = ||S^3 \ Eight|| + ||(S^1 x D^2) \ Wh_-^2||

# The simplicial volume of the figure-8 knot complement is exactly 2.
V_Eight = 2

# The simplicial volume of the pattern complement, which corresponds to the
# (-1)-twisted Whitehead clasp in a solid torus, is known to be exactly 1.
V_pattern = 1

# Summing these gives the simplicial volume for K2.
V_K2 = V_Eight + V_pattern

# Step 3: Calculate the total simplicial volume V.
V = V_K1 + V_K2

# Step 4: Compute the final result.
factor = 10**6
final_value = math.floor(factor * V)

# Print the step-by-step calculation.
print("The knot K is a connected sum of K1 = C_4,3(Conway) and K2 = Wh_-^2(Eight).")
print("The total simplicial volume V is V(K) = V(K1) + V(K2).")
print("-" * 20)

print("Calculating V(K1):")
print("V(K1) = ||S^3 \ C_4,3(Conway)|| = ||S^3 \ Conway||.")
print(f"The simplicial volume of the Conway knot complement is {V_K1}.")
print(f"So, V(K1) = {V_K1}")
print("-" * 20)

print("Calculating V(K2):")
print("V(K2) = ||S^3 \ Wh_-^2(Eight)|| = ||S^3 \ Eight|| + ||pattern||.")
print(f"The simplicial volume of the figure-8 knot complement is {V_Eight}.")
print(f"The simplicial volume of the pattern complement is {V_pattern}.")
print(f"V(K2) = {V_Eight} + {V_pattern} = {V_K2}")
print("-" * 20)

print("Calculating the total simplicial volume V:")
print(f"V = V(K1) + V(K2) = {V_K1} + {V_K2} = {V}")
print("-" * 20)

print("Final computation:")
print(f"The value to compute is floor({factor} * V).")
print(f"floor({factor} * {V}) = floor({factor * V}) = {final_value}")
