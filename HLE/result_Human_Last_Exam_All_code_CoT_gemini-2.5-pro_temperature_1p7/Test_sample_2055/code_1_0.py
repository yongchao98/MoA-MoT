import math

# Step 1: Define the distinct eigenvalues of the Seidel matrix S.
# These values are derived from the known eigenvalues of the Johnson scheme J(50,5).
# n = binomial(50,5) = 2118760
# k2 (degree of the graph) = 9900
# theta_0 = n - 1 - 2*k2 = 2118760 - 1 - 2*9900 = 2098959
# The other eigenvalues of S, theta_j, are -1 - 2*lambda_j, where lambda_j
# are the eigenvalues of the adjacency matrix A of the graph.
# The distinct eigenvalues of A (for j>0) are {3450, 450, -90, 30}.
# This leads to the following distinct eigenvalues for S (for j>0):
# -1 - 2*3450 = -6901
# -1 - 2*450 = -901
# -1 - 2*(-90) = 179
# -1 - 2*30 = -61
abs_theta_0 = 2098959
abs_theta_1 = 6901
abs_theta_2 = 901
abs_theta_3_4 = 179
abs_theta_5 = 61

# Step 2: The maximum order of an element in the Smith group is the largest
# invariant factor, d_n.
# Given that the prime factors of the absolute values of the distinct integer
# eigenvalues are disjoint, d_n is their product.

# List the numbers to be multiplied
values = [abs_theta_0, abs_theta_1, abs_theta_2, abs_theta_3_4, abs_theta_5]

# Step 3: Calculate the product.
result = 1
for v in values:
    result *= v

# Step 4: Print the final equation and the result.
equation_str = " * ".join(map(str, values))
print("The maximum order is the product of the absolute values of the distinct eigenvalues:")
print(f"{equation_str} = {result}")
