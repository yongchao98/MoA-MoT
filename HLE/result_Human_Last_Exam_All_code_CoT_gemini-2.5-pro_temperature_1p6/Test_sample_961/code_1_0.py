# Number of groups in the free product
num_groups = 19

# The power of each commutator
power = 30

# The stable commutator length of a single commutator [a, b] in a free group F(a,b)
scl_of_ci = 0.5

# Calculate the final stable commutator length
# scl(c) = sum_{i=1 to 19} scl(c_i^30) = sum_{i=1 to 19} 30 * scl(c_i)
#          = 19 * 30 * (1/2)
result = num_groups * power * scl_of_ci

# Print the equation with all the numbers
print(f"The stable commutator length is the result of the following calculation:")
print(f"{num_groups} * {power} * {scl_of_ci} = {result}")