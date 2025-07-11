import math
from scipy.special import wigner_3j

# Parameters identified from the source of the visualization.
# The nine plots correspond to j1 = 1, 2, ..., 9.
j2 = 8
j3 = 8
j1_values = range(1, 10)

infinity_norms = []

print("This script calculates the ratio of the maximum to the minimum infinity-norm")
print("for the nine Wigner 3-j symbols visualized in the image.")
print("-" * 50)
print("Calculating the infinity-norm (max absolute value) for each symbol...")

# Loop through each of the 9 symbols (matrices)
for j1 in j1_values:
    # This will store the maximum absolute value for the current j1
    current_max_abs_val = 0.0

    # Iterate over all possible m1 and m2 values.
    # The selection rules for m quantum numbers are: |m_i| <= j_i and m1+m2+m3=0.
    for m1 in range(-j1, j1 + 1):
        for m2 in range(-j2, j2 + 1):
            m3 = -m1 - m2

            # The wigner_3j function returns 0 if selection rules are not met.
            # We only need to ensure m3 is within its allowed range.
            if abs(m3) <= j3:
                val = wigner_3j(j1, j2, j3, m1, m2, m3)
                if abs(val) > current_max_abs_val:
                    current_max_abs_val = abs(val)

    infinity_norms.append(current_max_abs_val)
    print(f"  For j1 = {j1}, the norm is: {current_max_abs_val}")

# Find the overall maximum and minimum norms from the list
max_norm = max(infinity_norms)
min_norm = min(infinity_norms)

# Calculate the final ratio
# Check for min_norm being zero to avoid division by zero error
if min_norm == 0:
    ratio = float('inf')
else:
    ratio = max_norm / min_norm

print("-" * 50)
print(f"Overall maximum infinity-norm found: {max_norm}")
print(f"Overall minimum infinity-norm found: {min_norm}")
print("-" * 50)
print("Final ratio calculation:")
print(f"{max_norm} / {min_norm} = {ratio}")