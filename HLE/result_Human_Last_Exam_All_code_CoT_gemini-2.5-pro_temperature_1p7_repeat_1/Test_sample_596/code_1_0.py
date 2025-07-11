# Step 1: Define the order of the torsion subgroup.
# Based on the mathematical derivation, the first homology group of the Mobius band, H_1(M; Z), is Z (the integers).
# The torsion subgroup of Z is the trivial group {0}, whose order is 1.
order_of_torsion_subgroup = 1

# Step 2: Define the exponent from the theorem's formula.
exponent = 2

# Step 3: Calculate the final result according to the theorem |T(H_{d-1}(K))|^2.
result = order_of_torsion_subgroup ** exponent

# Step 4: Print the final equation with each number.
# The problem asks to output each number in the final equation.
# The equation is: (order_of_torsion_subgroup)^exponent = result.
# So we print the base, the exponent, and the final result.
print(f"The order of the torsion subgroup is: {order_of_torsion_subgroup}")
print(f"The exponent in the formula is: {exponent}")
print(f"The final number of non-collapsing rooted forests (in the weighted sense) is: {result}")
print(f"The final equation is: {order_of_torsion_subgroup}^{exponent} = {result}")
