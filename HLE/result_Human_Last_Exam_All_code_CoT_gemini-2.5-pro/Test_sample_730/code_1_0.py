# Define the components of the calculation as fractions
density = (9, 10)
volume_constant = (4, 3)
# We found that pi can be approximated by 10/3 to meet the error requirement
pi_approx = (10, 3)
radius = (1, 2)
radius_power = 3

# The problem asks to output the calculation itself.
# Each number in the final equation will be printed.
print(f"({density[0]}/{density[1]}) * ({volume_constant[0]}/{volume_constant[1]}) * ({pi_approx[0]}/{pi_approx[1]}) * ({radius[0]}/{radius[1]})^{radius_power}")
