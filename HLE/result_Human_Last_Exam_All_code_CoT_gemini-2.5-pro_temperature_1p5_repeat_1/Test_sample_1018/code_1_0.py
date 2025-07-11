import math

# Step 1: Define the genera and dimensions of the surfaces.
g1 = 31
g2 = 17
dim1 = 2
dim2 = 2

# Step 2: Calculate the simplicial volume of each surface using the formula ||S_g|| = 4*(g-1).
sv_g1 = 4 * (g1 - 1)
sv_g2 = 4 * (g2 - 1)

# Step 3: Calculate the binomial coefficient C(dim1+dim2, dim1) for the product formula.
# C(2+2, 2) = C(4, 2)
binomial_coeff = math.comb(dim1 + dim2, dim1)

# Step 4: Apply Gromov's product formula to get the final result.
# ||S_g1 x S_g2|| = C(4,2) * ||S_g1|| * ||S_g2||
final_volume = binomial_coeff * sv_g1 * sv_g2

# Print the final equation with all the computed values.
print(f"The simplicial volume is calculated using Gromov's product formula:")
print(f"||\u03A3_{g1} \u00D7 \u03A3_{g2}|| = C({dim1}+{dim2}, {dim1}) \u00D7 (4 \u00D7 ({g1}-1)) \u00D7 (4 \u00D7 ({g2}-1))")
print(f"Result = {binomial_coeff} \u00D7 {sv_g1} \u00D7 {sv_g2} = {final_volume}")