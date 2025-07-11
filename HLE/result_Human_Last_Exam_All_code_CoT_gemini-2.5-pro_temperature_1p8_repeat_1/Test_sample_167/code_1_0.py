# Define the parameters for the complete bipartite graph K_{m,n}
m = 1000
n = 1000

# The Alon-Tarsi number for a complete bipartite graph K_{m,n} is given by the formula min(m, n) + 1.
# We apply this formula to K_{1000,1000}.

# Calculate the minimum of m and n
min_val = min(m, n)

# Calculate the Alon-Tarsi number
alon_tarsi_number = min_val + 1

# Print the details of the calculation, showing each number in the equation.
print(f"To find the Alon-Tarsi number of K_{{{m},{n}}}, we use the formula AT(K_m,n) = min(m, n) + 1.")
print(f"Substituting the values m = {m} and n = {n}:")
print(f"AT(K_{{{m},{n}}}) = min({m}, {n}) + 1 = {min_val} + 1 = {alon_tarsi_number}")