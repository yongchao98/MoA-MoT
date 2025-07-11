import math

# Define the size of the generalized quaternion group.
group_size = 128

# The size of a generalized quaternion group is N = 2^n. We find n.
# n = log2(N)
n = int(math.log2(group_size))

# For a generalized quaternion group Q_{2^n} (with n >= 3),
# the exponent is 2^(n-1).
exponent = 2**(n - 1)
log_of_exponent = n - 1

# For a finite p-group G (where p is a prime), the number of power subgroups is
# given by the formula: 1 + log_p(exp(G)).
# In this case, the group is a 2-group, so p=2.
p = 2
num_power_subgroups = 1 + log_of_exponent

# Print the final result, showing each number in the calculation.
print(f"The number of power subgroups in the generalized quaternion group of size {group_size} is calculated as:")
print(f"1 + log{p}({exponent}) = 1 + {log_of_exponent} = {num_power_subgroups}")