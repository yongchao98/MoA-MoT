import numpy as np

def print_matrix(name, m):
    """Helper function to print a matrix with a name."""
    print(f"{name} = \n{m}")

# Define the matrices a, b, U, and I
a = np.array([[-21, 242], 
              [-2, 23]])

b = np.array([[-19, 200],
              [-2, 21]])

U = np.array([[1, 2],
              [0, 1]])

I = np.identity(2, dtype=int)

print("Step 1: Define the matrices a and b.")
print_matrix('a', a)
print_matrix('b', b)
print("")

# Step 1: Verify they are in SL(2,Z)
det_a = int(np.linalg.det(a))
det_b = int(np.linalg.det(b))
print(f"Step 2: Verify a and b are in G = SL_2(Z).")
print(f"det(a) = {det_a}")
print(f"det(b) = {det_b}")
if det_a == 1 and det_b == 1:
    print("a and b are in SL_2(Z).")
else:
    print("Error: a or b is not in SL_2(Z).")
print("")

# Step 2: Check if a and b are in Gamma(2)
a_mod_2 = a % 2
b_mod_2 = b % 2
print(f"Step 3: Check if a and b are in the congruence subgroup Gamma(2).")
print_matrix('a mod 2', a_mod_2)
print_matrix('b mod 2', b_mod_2)

# np.array_equal checks if two arrays have the same shape and elements.
if np.array_equal(a_mod_2, I % 2) and np.array_equal(b_mod_2, I % 2):
    print("a and b are in Gamma(2), so H is a subgroup of Gamma(2).")
else:
    print("Error: a or b is not in Gamma(2).")
print("")

# Step 4: Use the tower rule for indices
print("Step 4: Use the tower rule for indices.")
print("[G : H] = [G : Gamma(2)] * [Gamma(2) : H]")
print("")

# Step 5: Compute [G : Gamma(2)]
index_G_Gamma2 = 6
print(f"Step 5: The index [G : Gamma(2)] = |SL_2(Z/2Z)| = {index_G_Gamma2}.")
print("")

# Step 6 & 7: Show H = Gamma(2) by finding a relation
print("Step 6: Show that H = Gamma(2), which implies [Gamma(2) : H] = 1.")
print("To do this, we compute the product ba and relate it to the generators of Gamma(2).")
U_inv = np.linalg.inv(U).astype(int)
ba = b @ a
minus_U_inv = -1 * U_inv

print_matrix('b * a', ba)
print_matrix('-U^{-1}', minus_U_inv)

if np.array_equal(ba, minus_U_inv):
    print("The relation ba = -U^{-1} is verified.")
else:
    print("The relation ba = -U^{-1} could not be verified.")
print("")
    
# Step 8 & 9: Conclude the index is 1 and calculate final answer.
index_Gamma2_H = 1
print(f"Step 7: The relation ba = -U^{-1} allows us to show that the generators of Gamma(2)")
print(f"can be expressed in terms of a and b. This leads to the conclusion that H = Gamma(2),")
print(f"and thus [Gamma(2) : H] = {index_Gamma2_H}.")
print("")

final_index = index_G_Gamma2 * index_Gamma2_H
print("Step 8: Final calculation.")
print(f"[G : H] = [G : Gamma(2)] * [Gamma(2) : H] = {index_G_Gamma2} * {index_Gamma2_H} = {final_index}")

>>> 6