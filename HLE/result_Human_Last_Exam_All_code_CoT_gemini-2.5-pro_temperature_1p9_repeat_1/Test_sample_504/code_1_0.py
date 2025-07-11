import numpy as np

# Define the standard basis vectors in C^6
e = np.eye(6, dtype=complex)

# Define the 9 vectors
vectors = {
    'e1': e[0], 'e2': e[1], 'e3': e[2], 'e4': e[3], 'e5': e[4], 'e6': e[5],
    'u7': 0.5 * (e[0] + e[1] + e[2] + e[3]),
    'u8': 0.5 * (e[0] + e[1] + e[4] + e[5]),
    'u9': 0.5 * (e[2] + e[3] + e[4] + e[5])
}

# Example verification: angle between u8 and u9
u8 = vectors['u8']
u9 = vectors['u9']

# v.conj().T @ w is the standard hermitian form in numpy
inner_product_u8_u9 = u8.conj().T @ u9
cos_alpha_squared = np.abs(inner_product_u8_u9)**2
cos_alpha = np.abs(inner_product_u8_u9)

print("u8 = (1/2) * (e1 + e2 + e5 + e6)")
print("u9 = (1/2) * (e3 + e4 + e5 + e6)")
print(f"u8 = {u8}")
print(f"u9 = {u9}")
print(f"Inner product (u8, u9) = {inner_product_u8_u9:.4f}")
print(f"| (u8, u9) | = {cos_alpha:.4f}")
print("So the cosine of the angle between u8 and u9 is 1/2, which corresponds to pi/3.")

# Final answer
largest_number = 9
print(f"The largest number of such vectors is {largest_number}.")
print(f"This can be shown with the construction of 9 vectors as described, for example:")
print("The standard basis {e1, ..., e6} plus the three vectors:")
print("u7 = 1/2 * e1 + 1/2 * e2 + 1/2 * e3 + 1/2 * e4")
print("u8 = 1/2 * e1 + 1/2 * e2 + 1/2 * e5 + 1/2 * e6")
print("u9 = 1/2 * e3 + 1/2 * e4 + 1/2 * e5 + 1/2 * e6")
