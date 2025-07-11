import math

# Side lengths
s = 18

# --- Optimal Vertex Placement (conceptual) ---
# To avoid lattice points, we use small offsets epsilon and delta.
# For calculation, we use floating point approximations.
# A = (epsilon, delta)
# C = (s*sqrt(2) + epsilon, delta)
# B = (s/sqrt(2) + epsilon, s/sqrt(2) + delta)
# Let epsilon = 0.1, delta = 0.2

pA = (0.1, 0.2)
pB = (s * math.sqrt(2)/2 + 0.1, s * math.sqrt(2)/2 + 0.2)
pC = (s * math.sqrt(2) + 0.1, 0.2)

# --- Calculate squares for each side ---

# Side AC (A to C)
k_AC_nv = abs(math.floor(pC[0]) - math.floor(pA[0]))
k_AC_nh = abs(math.floor(pC[1]) - math.floor(pA[1]))
k_AC = 1 + k_AC_nv + k_AC_nh

# Side AB (A to B)
k_AB_nv = abs(math.floor(pB[0]) - math.floor(pA[0]))
k_AB_nh = abs(math.floor(pB[1]) - math.floor(pA[1]))
k_AB = 1 + k_AB_nv + k_AB_nh

# Side BC (B to C)
k_BC_nv = abs(math.floor(pC[0]) - math.floor(pB[0]))
k_BC_nh = abs(math.floor(pC[1]) - math.floor(pB[1]))
k_BC = 1 + k_BC_nv + k_BC_nh

# Total squares before removing overlaps
k_sum = k_AC + k_AB + k_BC

# Overlaps
# As determined by analysis of how segments with slope +/-1 and 0 interact near vertices.
# Overlap at A (AB and AC): 1 square
# Overlap at B (AB and BC): 1 square
# Overlap at C (BC and AC): 2 squares
overlaps = 1 + 1 + 2

# Final count
k_final = k_sum - overlaps

print(f"Number of squares for side AB: {k_AB}")
print(f"Number of squares for side BC: {k_BC}")
print(f"Number of squares for side AC: {k_AC}")
print(f"Sum of squares for all sides: {k_sum}")
print(f"Number of overlapping squares: {overlaps}")
print(f"The largest number of squares the perimeter can pass through is k:")
# The final equation requires each number to be output
# The calculation is k_sum - overlaps = k_final
print(f"{k_sum} - {overlaps} = {k_final}")
<<<73>>>