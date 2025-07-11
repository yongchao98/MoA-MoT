import numpy as np

# Step 1: Define the given 4-point DFTs from the problem.
# X(k) = [1, j, -1, -j]
# H(k) = [0, 1+j, 1, 1-j]
X = np.array([1, 1j, -1, -1j], dtype=np.complex128)
H = np.array([0, 1 + 1j, 1, 1 - 1j], dtype=np.complex128)

# The new sequence is an 8-point sequence y(n) formed by interleaving
# x(n) (even samples) and h(n) (odd samples). We need to find its DFT.
# The size of the resulting DFT is N=8.
N = 8

# Initialize the 8-point DFT Y(k) as an array of zeros.
Y = np.zeros(N, dtype=np.complex128)

# Step 2: Calculate the 8-point DFT Y(k) using the DIT-FFT synthesis formula.
# Y(k)     = X(k) + W_N^k * H(k)
# Y(k+N/2) = X(k) - W_N^k * H(k)
# We compute this for k = 0, 1, 2, 3.

# Calculate the twiddle factors W_8^k for k = 0, 1, 2, 3.
k_half = np.arange(N // 2)
W = np.exp(-2j * np.pi * k_half / N)

# Calculate the term T(k) = W_8^k * H(k)
T = W * H

# Apply the butterfly formulas
Y[:N//2] = X + T
Y[N//2:] = X - T

# Step 3: The problem asks for the "4-point DFT" of the 8-point sequence.
# We interpret this as decimating the 8-point DFT Y(k) in frequency,
# by taking every second sample to produce a 4-point sequence Z(m).
# Z(m) = Y(2m) for m = 0, 1, 2, 3.
Z = Y[::2]

# Step 4: Print the final 4-point DFT Z(m).
# The format fulfills the requirement "output each number in the final equation".
print("The resulting 4-point DFT, Z(m), is:")
for m in range(len(Z)):
    # The calculated values have integer real and imaginary parts.
    val = Z[m]
    print(f"Z({m}) = {val}")
