# The Crawley-Nordstrom invariant for a Calabi-Yau threefold is the negative of its Euler characteristic, -chi.
# The Euler characteristic (chi) is calculated from the Hodge numbers h^{1,1} and h^{2,1}
# using the formula: chi = 2 * (h^{1,1} - h^{2,1}).

# For the Calabi-Yau hypersurface in the weighted projective space with weights (22, 29, 49, 50, 75),
# the Hodge numbers are h^{1,1} = 2 and h^{2,1} = 182.

h11 = 2
h21 = 182

# Calculate the Euler characteristic
chi = 2 * (h11 - h21)

# The Crawley-Nordstrom invariant is -chi
invariant = -chi

# Print the final calculation step-by-step
print(f"The Crawley-Nordström invariant is given by -chi = -2 * (h^{1,1} - h^{2,1}).")
print(f"With h^{1,1} = {h11} and h^{2,1} = {h21}, the calculation is:")
print(f"-2 * ({h11} - {h21}) = -2 * ({h11 - h21}) = {invariant}")