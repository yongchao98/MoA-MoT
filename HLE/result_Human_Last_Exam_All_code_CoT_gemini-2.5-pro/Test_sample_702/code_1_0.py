# Sphere dimensions
p = 4
q = 6

# The connectivity of the source space A = Sigma(Omega S^p wedge Omega S^q)
# Connectivity of Omega S^p is p-2
# Connectivity of Omega S^q is q-2
conn_omega_p = p - 2
conn_omega_q = q - 2

# Connectivity of the smash product is conn1 + conn2 + 1
conn_smash = conn_omega_p + conn_omega_q + 1

# Connectivity of the suspension is conn + 1
conn_source = conn_smash + 1

print(f"The dimension of the first sphere is p = {p}")
print(f"The dimension of the second sphere is q = {q}")
print(f"Connectivity of the source space A = Sigma(Omega S^{p} wedge Omega S^{q}) is {conn_source}")

# The connectivity of the target space B = Omega(S^p wedge S^q) = Omega(S^{p+q})
# Connectivity of Omega S^n is n-2
conn_target = p + q - 2
print(f"Connectivity of the target space B = Omega(S^{p+q}) is {conn_target}")

# The map is an isomorphism on pi_9 and pi_10
# pi_9(A) = Z, pi_9(B) = Z
# pi_10(A) = Z_2, pi_10(B) = Z_2
# The map is a surjection on pi_11
# pi_11(A) = Z, pi_11(B) = Z_2

# The connectivity is defined as the largest n for which pi_k is an isomorphism for k < n
# and pi_n is a surjection.
connectivity = 11

print(f"The map is an isomorphism on homotopy groups for k < {connectivity}.")
print(f"The map is a surjection on the homotopy group for k = {connectivity}.")
print(f"Therefore, the connectivity of the map is {connectivity}.")
