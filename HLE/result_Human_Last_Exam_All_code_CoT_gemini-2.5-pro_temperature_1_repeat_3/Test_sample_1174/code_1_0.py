# The exponents determined from the physical analysis of the system.
n1 = 1
n2 = 1
n3 = -2
n4 = 0
n5 = -2
n6 = -1

# The coefficients k from the summation.
k1, k2, k3, k4, k5, k6 = 1, 2, 3, 4, 5, 6

# Calculate the sum.
total_sum = k1*n1 + k2*n2 + k3*n3 + k4*n4 + k5*n5 + k6*n6

# Construct the equation string.
equation_str = (
    f"{k1}*({n1}) + {k2}*({n2}) + {k3}*({n3}) + "
    f"{k4}*({n4}) + {k5}*({n5}) + {k6}*({n6}) = {total_sum}"
)

# Print the final result in the requested format.
print("The final equation is:")
print(equation_str)
print("\nThe value of the sum is:")
print(total_sum)