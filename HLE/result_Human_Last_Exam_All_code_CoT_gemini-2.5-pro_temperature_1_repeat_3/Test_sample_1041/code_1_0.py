import math

# Define the parameters of the problem
R = 1000.0  # Radius of the disk
S = (0, 300)   # Starting point
T1 = (0, 0)    # First target point
T2 = (2, 0)    # Second target point
a = 1.0        # Lattice spacing, used as the effective radius of a single point target

# Calculate the required distances
dist_S_T1 = math.sqrt((S[0] - T1[0])**2 + (S[1] - T1[1])**2)
dist_S_T2 = math.sqrt((S[0] - T2[0])**2 + (S[1] - T2[1])**2)
dist_T1_T2 = math.sqrt((T1[0] - T2[0])**2 + (T1[1] - T2[1])**2)

# The terms in the formula correspond to Green's function approximations (ignoring a constant factor)
# Numerator terms
g_S_T1 = math.log(R / dist_S_T1)
g_S_T2 = math.log(R / dist_S_T2)

# Denominator terms
g_T1_T1 = math.log(R / a)
g_T1_T2 = math.log(R / dist_T1_T2)

# Calculate the final probability
# P = (g_S_T1 + g_S_T2) / (g_T1_T1 + g_T1_T2)
numerator = g_S_T1 + g_S_T2
denominator = g_T1_T1 + g_T1_T2
probability = numerator / denominator

# Output the components of the calculation
print("This script calculates the probability based on the continuous approximation of the random walk.")
print("The formula used is: P = [ln(R/|S-T1|) + ln(R/|S-T2|)] / [ln(R/a) + ln(R/|T1-T2|)]\n")

print("--- Parameters ---")
print(f"Disk Radius (R): {R}")
print(f"Starting Point (S): {S}")
print(f"Target Set (T): {{{T1}, {T2}}}")
print(f"Lattice Spacing (a): {a}\n")

print("--- Calculated Distances ---")
print(f"|S - T1|: {dist_S_T1}")
print(f"|S - T2|: {dist_S_T2:.6f}")
print(f"|T1 - T2|: {dist_T1_T2}\n")

print("--- Numerator and Denominator of the Formula ---")
print(f"Numerator = ln({R}/{dist_S_T1}) + ln({R}/{dist_S_T2:.6f}) = {g_S_T1:.6f} + {g_S_T2:.6f} = {numerator:.6f}")
print(f"Denominator = ln({R}/{a}) + ln({R}/{dist_T1_T2}) = {g_T1_T1:.6f} + {g_T1_T2:.6f} = {denominator:.6f}\n")

print("--- Final Probability ---")
print(f"P = {numerator:.6f} / {denominator:.6f} = {probability:.6f}")
print(f"The probability rounded to three significant digits is: {probability:.3g}")
