# The closed-form expression for the coefficient a_{2n+1} for n>=1 is 0.
expression_a_2n_plus_1 = "0"

# The closed-form expression for the coefficient a_{2n} for n>=1 is 2^(2n-1) * ((n-1)!)^2 / (2n)!.
# We represent it here using Python-like syntax for clarity.
expression_a_2n = "2**(2*n - 1) * ((n - 1)!)**2 / (2*n)!"

# Print the two expressions separated by a comma, as requested.
print(f"{expression_a_2n_plus_1}, {expression_a_2n}")