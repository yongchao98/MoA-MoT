import numpy as np

m = 4

# Our successful construction of n = m+1 vectors
V = np.identity(m)
all_ones = np.ones((m, 1))
V = np.hstack([V, all_ones])

print(f"--- Checking our set of n = m+1 = {m+1} vectors for m={m} ---")
print("This set is {e_1, e_2, ..., e_m, all_ones_vector}")
print("It is in general linear position. For instance, the determinant of {e_1, e_2, e_3, all_ones_vector} is:")
# Example subset of size m
subset_1 = V[:, [0, 1, 2, 4]]
det_1 = np.linalg.det(subset_1)
print(det_1)
print("-" * 20)


# Now, let's try to add an (m+2)th vector
# v_{m+2} = e_1 + e_2 = (1,1,0,0)
v_new = np.zeros((m, 1))
v_new[0, 0] = 1
v_new[1, 0] = 1

print(f"--- Trying to add a {m+2}-th vector: v = e_1 + e_2 ---")
print("Consider the subset of m vectors {e_1, e_2, e_3, v}")

# This subset is {e_1, e_2, e_3, e_1+e_2} for m=4
subset_2_vectors = {
    'e1': V[:, 0],
    'e2': V[:, 1],
    'e3': V[:, 2],
    'v_new': v_new.flatten()
}

subset_2_matrix = np.array([V[:,0], V[:,1], V[:,2], v_new.flatten()]).T
print("The matrix formed by these vectors is:")
print(subset_2_matrix)

det_2 = np.linalg.det(subset_2_matrix)
print(f"\nThe determinant of this matrix is: {det_2}")
if np.isclose(det_2, 0):
    print("This subset is linearly dependent, which violates the condition.")
else:
    print("This subset is linearly independent.")
    
print("\nWe can write the dependency explicitly:")
c1 = 1
c2 = 1
c3 = -1
v1 = subset_2_vectors['e1']
v2 = subset_2_vectors['e2']
v_new_vec = subset_2_vectors['v_new']

result_vector = c1 * v1 + c2 * v2 + c3 * v_new_vec

# Print the final equation with each number
print(f"{c1} * {tuple(v1)} + {c2} * {tuple(v2)} + {c3} * {tuple(v_new_vec)} = {tuple(result_vector)}")
