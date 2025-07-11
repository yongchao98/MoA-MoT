import numpy as np

# Suppose that X(k), H(k) are 4-point DFT for 4-point sequences x(n), h(n), respectively.
# X(k)= [1, j, -1, -j].
# H(k)=[0,1+j,1,1-j].
# Find 4-point DFT for sequence {x(0),h(0),x(1),h(1),x(2),h(2),x(3),h(3)}

# --- Input Data ---
# Given 4-point DFTs
X = np.array([1, 1j, -1, -1j])
H = np.array([0, 1 + 1j, 1, 1 - 1j])

# --- Theory ---
# The new 8-point sequence is y(n) = {x(0),h(0),x(1),h(1),x(2),h(2),x(3),h(3)}.
# The 4-point DFT of y(n), let's call it Y_4(k), is found by sampling the 8-point DFT of y(n), Y_8(k).
# The relation is: Y_4(k) = Y_8(2k) for k = 0, 1, 2, 3.
# The 8-point DFT Y_8(k) can be expressed in terms of X(k) and H(k) as:
# Y_8(k) = X(k mod 4) + W_8^k * H(k mod 4), where W_8 = exp(-2j*pi/8).

# We need to compute:
# Y_4(0) = Y_8(0) = X(0) + W_8^0 * H(0)
# Y_4(1) = Y_8(2) = X(2) + W_8^2 * H(2)
# Y_4(2) = Y_8(4) = X(0) + W_8^4 * H(0)
# Y_4(3) = Y_8(6) = X(2) + W_8^6 * H(2)

# --- Calculation ---
# Calculate the required twiddle factors W_8^n
W_8_0 = np.exp(-2j * np.pi * 0 / 8)  # Equals 1
W_8_2 = np.exp(-2j * np.pi * 2 / 8)  # Equals -j
W_8_4 = np.exp(-2j * np.pi * 4 / 8)  # Equals -1
W_8_6 = np.exp(-2j * np.pi * 6 / 8)  # Equals j

# Calculate the 4 points of the resulting DFT, Y_4(k)
Y4_0 = X[0] + W_8_0 * H[0]
Y4_1 = X[2] + W_8_2 * H[2]
Y4_2 = X[0] + W_8_4 * H[0]
Y4_3 = X[2] + W_8_6 * H[2]

Y4 = np.array([Y4_0, Y4_1, Y4_2, Y4_3])

# --- Output ---
print("Calculating the 4-point DFT, Y_4(k):")
print("-" * 40)
print(f"Y_4(0) = X(0) + (W_8^0) * H(0)")
print(f"       = {X[0]} + ({np.round(W_8_0, 4)}) * ({H[0]})")
print(f"       = {np.round(Y4_0, 4)}\n")

print(f"Y_4(1) = X(2) + (W_8^2) * H(2)")
print(f"       = ({X[2]}) + ({np.round(W_8_2, 4)}) * ({H[2]})")
print(f"       = {np.round(Y4_1, 4)}\n")

print(f"Y_4(2) = X(0) + (W_8^4) * H(0)")
print(f"       = {X[0]} + ({np.round(W_8_4, 4)}) * ({H[0]})")
print(f"       = {np.round(Y4_2, 4)}\n")

print(f"Y_4(3) = X(2) + (W_8^6) * H(2)")
print(f"       = ({X[2]}) + ({np.round(W_8_6, 4)}) * ({H[2]})")
print(f"       = {np.round(Y4_3, 4)}\n")

print("-" * 40)
print("The final 4-point DFT is:")
# Use a custom formatting to avoid parentheses for simple numbers
formatted_result = [f"{v.real if v.imag == 0 else v}" for v in np.round(Y4, 4)]
print(f"[{', '.join(formatted_result)}]")
