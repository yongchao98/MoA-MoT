import numpy as np

# Observed intensities from the table for the molecular ion cluster
I_M = 29.0   # Intensity at m/z 262
I_M_plus_2 = 12.9 # Intensity at m/z 264
I_M_plus_4 = 3.22 # Intensity at m/z 266

# Theoretical ratio for a compound with two Cl atoms is approximately 9:6:1
# We can normalize the observed intensities to the M+4 peak
# and see how well they match the theoretical 9:6:1 ratio.
if I_M_plus_4 > 0:
    observed_ratio_M = I_M / I_M_plus_4
    observed_ratio_M_plus_2 = I_M_plus_2 / I_M_plus_4
else:
    observed_ratio_M = float('inf')
    observed_ratio_M_plus_2 = float('inf')

theoretical_M = 9.0
theoretical_M_plus_2 = 6.0
theoretical_M_plus_4 = 1.0

print("--- Isotope Pattern Analysis for a Dichloro Compound ---")
print(f"Observed Peaks: m/z 262 (I={I_M}), m/z 264 (I={I_M_plus_2}), m/z 266 (I={I_M_plus_4})")
print(f"Theoretical Ratio (M:M+2:M+4) for Cl2: {theoretical_M}:{theoretical_M_plus_2}:{theoretical_M_plus_4}")
print(f"Observed Ratio (normalized to M+4): {observed_ratio_M:.2f}:{observed_ratio_M_plus_2:.2f}:1.0")

# Also check the M-Cl fragment at m/z 227 and its isotope at m/z 229
I_227 = 66.7
I_229 = 22.2
# Theoretical ratio for one Cl atom is 3:1
observed_ratio_227_229 = I_227 / I_229 if I_229 > 0 else float('inf')
print("\n--- Analysis of Fragment with One Chlorine Atom ---")
print(f"Fragment Peaks: m/z 227 (I={I_227}), m/z 229 (I={I_229})")
print(f"Theoretical Ratio for a single Cl atom: 3:1")
print(f"Observed Ratio (227:229): {observed_ratio_227_229:.2f}:1")
