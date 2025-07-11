# Number of free groups in the free product
num_groups = 19

# The power to which each commutator is raised
power = 30

# The stable commutator length of a single commutator [a,b] in a free group F_2 is 1.
# The scl of c_i^30 in F_i is power * scl(c_i) = 30 * 1 = 30.
scl_per_term = power

# The total scl is the sum of the scl of each term, as per Burger's theorem.
total_scl = num_groups * scl_per_term

# Output the numbers in the final equation.
print(f"The final calculation is the product of the number of terms and the scl of each term.")
print(f"{num_groups} * {scl_per_term} = {total_scl}")