# Define the dimensions of the spheres
n1 = 4
n2 = 6

# Step 1: Calculate the connectivity of the loop spaces
# The connectivity of Omega S^n is n-2 for n >= 2
conn_Omega_S_n1 = n1 - 2
conn_Omega_S_n2 = n2 - 2
print(f"Connectivity of Omega S^{n1}: {conn_Omega_S_n1}")
print(f"Connectivity of Omega S^{n2}: {conn_Omega_S_n2}")

# Step 2: Calculate the connectivity of the smash product of the loop spaces
# The connectivity of A wedge B is conn(A) + conn(B) + 1
conn_smash = conn_Omega_S_n1 + conn_Omega_S_n2 + 1
print(f"Connectivity of (Omega S^{n1} wedge Omega S^{n2}): {conn_Omega_S_n1} + {conn_Omega_S_n2} + 1 = {conn_smash}")

# Step 3: Calculate the connectivity of the source space (suspension of the smash product)
# The connectivity of Sigma X is conn(X) + 1
conn_source = conn_smash + 1
print(f"Connectivity of source space Sigma(Omega S^{n1} wedge Omega S^{n2}): {conn_smash} + 1 = {conn_source}")

# Step 4: Calculate the connectivity of the target space
# The target is Omega(S^n1 wedge S^n2) = Omega(S^(n1+n2))
n_target = n1 + n2
conn_target = n_target - 2
print(f"Connectivity of target space Omega(S^{n1+n2}): {n_target} - 2 = {conn_target}")

# Step 5: Calculate the connectivity of the fiber of the null-homotopic map
# The map is null-homotopic. Its fiber is X x Omega(Y).
# X is the source space, which is 8-connected.
# Y is the target space, Omega(S^10).
# We need the connectivity of Omega(Y) = Omega(Omega(S^10)) = Omega^2(S^10).
# Connectivity of Omega^2(S^k) is k-3.
conn_Omega_Y = n_target - 3
print(f"Connectivity of Omega(Target) = Omega^2(S^{n_target}): {n_target} - 3 = {conn_Omega_Y}")

# The connectivity of the fiber is min(conn(Source), conn(Omega(Target)))
conn_fiber = min(conn_source, conn_Omega_Y)
print(f"Connectivity of the fiber: min({conn_source}, {conn_Omega_Y}) = {conn_fiber}")

# Step 6: The connectivity of the map is conn(fiber) + 1
conn_map = conn_fiber + 1
print(f"The connectivity of the map is conn(fiber) + 1 = {conn_fiber} + 1 = {conn_map}")
