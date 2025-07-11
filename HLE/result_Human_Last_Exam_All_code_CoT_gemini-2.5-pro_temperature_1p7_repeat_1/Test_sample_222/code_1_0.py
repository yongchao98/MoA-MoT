import numpy as np

# This script calculates the 8-point DFT of an interleaved sequence y(n)
# given the 4-point DFTs of the original sequences, X(k) and H(k).
#
# The interleaved sequence is y(n) = {x(0), h(0), x(1), h(1), ...}
# The formula derived from DFT properties (upsampling, time-shifting, linearity) is:
# Y(k) = X(k mod 4) + (W_8^k) * H(k mod 4)
# where W_8 = exp(-j*2*pi/8) is the twiddle factor.

# Step 1: Define the given 4-point DFTs as numpy arrays.
X_k = np.array([1, 1j, -1, -1j])
H_k = np.array([0, 1 + 1j, 1, 1 - 1j])

print("Given 4-point DFTs:")
print(f"X(k) = {X_k}")
print(f"H(k) = {H_k}\n")

# Step 2: Set up parameters for the calculation.
N = 4  # Size of original DFTs
M = 8  # Size of the target DFT

# According to the upsampling property, the 8-point DFTs of the upsampled
# sequences are periodic repetitions of the original 4-point DFTs.
G1_k = np.tile(X_k, M // N)  # This represents X(k mod 4) for k=0..7
G2_k = np.tile(H_k, M // N)  # This represents H(k mod 4) for k=0..7

# Calculate the twiddle factors W_M^k needed for the time-shift property.
k_values = np.arange(M)
W_M_k = np.exp(-2j * np.pi * k_values / M)

# Step 3: Calculate and print the result for each value of k from 0 to 7.
print("Calculating the 8-point DFT Y(k) using the formula:")
print(f"Y(k) = X(k mod {N}) + (W_{M}^k) * H(k mod {N})")
print("=" * 60)

final_Y_k = []
for k in range(M):
    # Get the components for the current value of k
    g1 = G1_k[k]
    g2 = G2_k[k]
    w = W_M_k[k]

    # Calculate the final result for Y(k)
    y = g1 + w * g2
    final_Y_k.append(y)
    
    # Format the numbers for a clean, readable output equation
    g1_str = f"({g1:.2f})".replace(".00", "").replace("+0.00j", "").replace("-0.00j", "").replace("0.00-", "-").replace("0.00+", "")
    w_str = f"({w.real:.3f}{w.imag:+.3f}j)"
    g2_str = f"({g2:.2f})".replace(".00", "").replace("+0.00j", "").replace("-0.00j", "").replace("0.00-", "-").replace("0.00+", "")
    y_str = f"({y.real:.4f}{y.imag:+.4f}j)"

    # Print the equation with all the numbers for the current step
    print(f"Y({k}) = {g1_str} + {w_str} * {g2_str}")
    print(f"     = {y_str}\n")

# Convert list to a numpy array for final display
Y_k_result = np.array(final_Y_k)

print("=" * 60)
print("The final 8-point DFT sequence Y(k) is:")
# We use np.round to present a clean final result
print(np.round(Y_k_result, 4))
<<<[ 1.    +0.j      1.4142+1.j     -1.    -1.j     -1.4142-1.j
  1.    +0.j     -1.4142+1.j     -1.    +1.j      1.4142-1.j    ]>>>