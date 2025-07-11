# Set the dimensions of the spheres
p = 4
q = 6

# The connectivity of the map Sigma(Omega S^p wedge Omega S^q) -> Omega(S^p wedge S^q)
# is given by the formula: p + q + min(p, q) - 3.
# S^p wedge S^q is homotopy equivalent to S^{p+q}.

print("The connectivity of the map is calculated using the formula: p + q + min(p, q) - 3.")
print(f"Here, p = {p} and q = {q}.")
print("")

# Calculate the minimum of p and q
min_val = min(p, q)

# Calculate the connectivity using the formula
connectivity = p + q + min_val - 3

# Print the step-by-step calculation
print("The final equation is:")
print(f"Connectivity = {p} + {q} + min({p}, {q}) - 3")
print(f"Connectivity = {p} + {q} + {min_val} - 3")
print(f"Connectivity = {p + q + min_val} - 3")
print(f"Connectivity = {connectivity}")

# The final answer in the required format is derived from this calculation.
# print(f"\n<<<{connectivity}>>>")