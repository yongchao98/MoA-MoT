# Define the dimensions of the spheres in the problem.
p = 4
q = 6

print(f"We want to find the connectivity of the map f: Sigma(Omega S^{p} wedge Omega S^{q}) -> Omega(S^{p} wedge S^{q}).")
print("This map fits into a homotopy fibration sequence where its fiber is Omega(S^p * S^q).")
print("The connectivity of the map is determined by the connectivity of its fiber.")
print("-" * 20)

# Step 1: Calculate the dimension of the join of the two spheres.
# The join of S^p and S^q is S^(p+q+1).
dim_join = p + q + 1
print(f"The join of S^{p} and S^{q} is S^({p}+{q}+1) = S^{dim_join}.")

# Step 2: Calculate the connectivity of the join sphere.
# A sphere S^n is (n-1)-connected.
conn_join_sphere = dim_join - 1
print(f"The sphere S^{dim_join} is ({dim_join}-1) = {conn_join_sphere}-connected.")

# Step 3: Calculate the connectivity of the fiber of the map.
# The fiber is the loop space of the join sphere, Omega(S^dim_join).
# The connectivity of Omega(X) is one less than the connectivity of X.
conn_fiber = conn_join_sphere - 1
print(f"The fiber of the map, Omega(S^{dim_join}), is therefore ({conn_join_sphere}-1) = {conn_fiber}-connected.")

# Step 4: Calculate the connectivity of the map itself.
# A map is k-connected if its fiber is (k-1)-connected.
# So, map_connectivity = fiber_connectivity + 1.
conn_map = conn_fiber + 1
print("-" * 20)
print(f"The connectivity of the map is the connectivity of its fiber plus one.")
print(f"Therefore, the final connectivity is {conn_fiber} + 1 = {conn_map}.")
