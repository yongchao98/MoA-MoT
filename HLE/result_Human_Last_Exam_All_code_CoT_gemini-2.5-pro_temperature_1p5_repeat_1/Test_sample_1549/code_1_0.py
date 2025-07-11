# The space is X = [0,1]^3, a product of 3 spaces.
k = 3

# The compactness number for each factor space X_i = [0,1] is 2.
# We create a list to hold the compactness number for each of the k factor spaces.
n_values = [2] * k

# The formula for the compactness number of a product space is:
# [Product(X_i)] = Sum([X_i] - 1) for i=1 to k, all plus 1.
# We calculate each ([X_i] - 1) term.
terms = [n - 1 for n in n_values]

# Sum the terms and add 1 to get the final result.
result = sum(terms) + 1

# We want to print the full equation for clarity.
# Create the string representation of the sum.
equation_parts = [f"({n}-1)" for n in n_values]
equation_str = " + ".join(equation_parts) + " + 1"

print(f"To find [X] for X = [0,1]^3, we use the product formula.")
print(f"Given that for a single interval I = [0,1], we have [I] = 2.")
print(f"The formula for a product of k spaces is Sum([X_i] - 1) + 1.")
print(f"For [0,1]^3, the calculation is:")
print(f"{equation_str} = {result}")
