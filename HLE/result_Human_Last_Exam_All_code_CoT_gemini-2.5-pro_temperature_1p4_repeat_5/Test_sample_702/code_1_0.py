# The dimensions of the spheres in the problem definition
m = 4
n = 6

# The connectivity of the source space, Sigma(Omega S^m wedge Omega S^n)
# Connectivity of Omega S^m is m-2
# Connectivity of Omega S^n is n-2
# Connectivity of the smash product is (m-2) + (n-2) + 1 = m+n-3
# Connectivity of the suspension is (m+n-3)+1 = m+n-2
connectivity_source = m + n - 2
print(f"Let m={m} and n={n}.")
print(f"The source space Sigma(Omega S^{m} wedge Omega S^{n}) is {connectivity_source}-connected.")

# The connectivity of the target space, Omega(S^m wedge S^n) = Omega(S^{m+n})
# Connectivity is (m+n)-2
connectivity_target = m + n - 2
print(f"The target space Omega(S^{m} wedge S^{n}) is {connectivity_target}-connected.")

# The map is an isomorphism on pi_k for k <= m+n-1
iso_up_to = m + n - 1
print(f"The map is an isomorphism on homotopy groups up to dimension {iso_up_to}.")

# The connectivity of the map is the largest n such that pi_k is an
# isomorphism for k < n and an epimorphism for k = n.
# Since the map is an isomorphism for k <= m+n-1 = 9, it is at least 10-connected.
# At dimension m+n=10, the map is from Z_2 to Z_2, which is not an epimorphism.
# Thus, the connectivity is m+n.
connectivity_map = m + n
print(f"The map is an epimorphism on pi_{m+n-1}, but not on pi_{m+n}.")
print(f"The connectivity of the map is m + n = {m} + {n} = {connectivity_map}.")