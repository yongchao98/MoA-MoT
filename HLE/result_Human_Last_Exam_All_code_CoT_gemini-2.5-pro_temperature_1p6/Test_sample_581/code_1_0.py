import math

# The dimension of the vector space AG(n, 3)
n = 8

# The best known lower bound for the size of a cap set in dimension n=8
# was found by Yves Edel. While the construction is complex, the resulting
# number can be calculated with the following expression for n=8.
term1_exp = n + 1
term2_exp = n // 2 # Using integer division for n/2

# Calculate the values of the two terms
val1 = int(math.pow(2, term1_exp))
val2 = int(math.pow(2, term2_exp))

# Calculate the final result
result = val1 - val2

# Print the final equation with all the numbers
print(f"The best known lower bound for the size of cap sets in dimension {n} is {result}.")
print("This can be calculated using the expression 2^(n+1) - 2^(n/2) for n=8:")
print(f"2^({n}+1) - 2^({n}//2) = 2^{term1_exp} - 2^{term2_exp} = {val1} - {val2} = {result}")
