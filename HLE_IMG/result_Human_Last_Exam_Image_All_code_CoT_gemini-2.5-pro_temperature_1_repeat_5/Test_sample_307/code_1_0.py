import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import N

def calculate_infinity_norm(j1, j2, j3):
    """
    Calculates the infinity-norm for a Wigner 3-j symbol matrix.
    The norm is the maximum absolute row sum.
    """
    m2_range = np.arange(-j2, j2 + 1)
    m1_range = np.arange(-j1, j1 + 1)
    
    max_row_sum = 0
    
    # Iterate over rows (m2 values)
    for m2 in m2_range:
        current_row_sum = 0
        # Sum up absolute values across the columns (m1 values)
        for m1 in m1_range:
            m3 = -m1 - m2
            # The wigner_3j function returns 0 if any selection rule is violated.
            # We convert the symbolic result from sympy to a float.
            val = float(N(wigner_3j(j1, j2, j3, m1, m2, m3)))
            current_row_sum += abs(val)
        
        if current_row_sum > max_row_sum:
            max_row_sum = current_row_sum
            
    return max_row_sum

# Based on the analysis of the image caption and fixing the inconsistencies,
# these are the most likely (j1, j2, j3) parameters for the nine plots.
params = [
    (4, 2, 2),  # Plot 1: (8,2,2) is invalid, (6,2,2) is invalid, (4,2,2) is valid.
    (6, 2, 4),  # Plot 2: (8,2,4) is invalid, (6,2,4) is valid.
    (8, 2, 6),  # Plot 3: (8,2,6) is valid.
    (6, 4, 2),  # Plot 4: (8,4,2) is invalid, (6,4,2) is valid.
    (8, 4, 4),  # Plot 5: (8,4,4) is valid.
    (8, 4, 6),  # Plot 6: (8,4,6) is valid.
    (8, 6, 2),  # Plot 7: (8,6,2) is valid.
    (8, 6, 4),  # Plot 8: (8,6,4) is valid.
    (8, 6, 6)   # Plot 9: (8,6,6) is valid.
]

norms = []
print("Calculating infinity-norm for each Wigner 3-j symbol matrix...")
for i, (j1, j2, j3) in enumerate(params):
    norm = calculate_infinity_norm(j1, j2, j3)
    norms.append(norm)
    print(f"Plot {i+1} with (j1,j2,j3) = {j1,j2,j3}: Infinity-norm = {norm:.6f}")

# Find the maximum and minimum norms
max_norm = max(norms)
min_norm = min(norms)

print("\n--------------------------------------------------")
print(f"Maximum norm found: {max_norm}")
print(f"Minimum norm found: {min_norm}")

# Calculate the final ratio
ratio = max_norm / min_norm

# Output the final equation with all the numbers
print("\nFinal Ratio Calculation:")
print(f"{max_norm} / {min_norm} = {ratio}")
print("--------------------------------------------------")
