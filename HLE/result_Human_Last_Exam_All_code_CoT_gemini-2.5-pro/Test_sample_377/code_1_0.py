import math

# Define the parameters for the extraspecial group S = 3^(1+2)_+
p = 3
n = 1

# The number of blocks of kG is equal to the number of conjugacy classes of S.
# The number of conjugacy classes of an extraspecial group of order p^(1+2n)
# is given by the formula: p^(2n) + p - 1.

# Calculate the number of classes.
num_classes_S = p**(2 * n) + p - 1

# Print the reasoning and the final calculation.
print("The group G is 2-solvable, so the number of blocks equals the number of conjugacy classes of 2'-elements.")
print("These classes are in one-to-one correspondence with the conjugacy classes of the subgroup S.")
print("S is an extraspecial group of order 27, i.e., p^(1+2n) with p=3 and n=1.")
print("The number of conjugacy classes is calculated using the formula p^(2n) + p - 1.")
print(f"Calculation: {p}^(2*{n}) + {p} - 1 = {int(p**(2*n))} + {p} - 1 = {int(num_classes_S)}")
print(f"Therefore, the number of blocks of kG is {int(num_classes_S)}.")
