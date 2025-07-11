# Values of the exponents based on physical analysis
n1 = 1
n2 = 1
n3 = -1
n4 = 0
n5 = -2
n6 = -1.5

# The expression to be calculated is sum_{k=1 to 6} k * n_k
# Let's compute this sum.
# The variable names k1, k2, ... correspond to the weights in the sum
k1, k2, k3, k4, k5, k6 = 1, 2, 3, 4, 5, 6

term1 = k1 * n1
term2 = k2 * n2
term3 = k3 * n3
term4 = k4 * n4
term5 = k5 * n5
term6 = k6 * n6

total_sum = term1 + term2 + term3 + term4 + term5 + term6

# We print the full equation as requested, showing each number.
print(f"The exponents are: n1={n1}, n2={n2}, n3={n3}, n4={n4}, n5={n5}, n6={n6}")
print("The calculation is:")
print(f"({k1} * {n1}) + ({k2} * {n2}) + ({k3} * {n3}) + ({k4} * {n4}) + ({k5} * {n5}) + ({k6} * {n6}) = {total_sum}")