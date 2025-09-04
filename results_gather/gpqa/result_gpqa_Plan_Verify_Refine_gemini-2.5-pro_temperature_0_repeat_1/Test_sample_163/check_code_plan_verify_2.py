import math

# Given data from the problem
# System 1
P1 = 2.0  # years
K1a = 10.0 # km/s
K1b = 5.0  # km/s

# System 2
P2 = 1.0  # year
K2a = 15.0 # km/s
K2b = 10.0 # km/s

# The problem states the systems are eclipsing, which implies the inclination angle i is ~90 degrees.
# For i = 90 degrees, sin(i) = 1.
# The formula for the total mass of a binary system is derived from Kepler's Third Law:
# M_total = (P * (K1 + K2)**3) / (2 * pi * G * (sin(i))**3)
# Since we are calculating a ratio of masses, the constant (2 * pi * G) cancels out.
# We assume sin(i) = 1 for both systems as they are eclipsing.
# Ratio = M1_total / M2_total = (P1 * (K1a + K1b)**3) / (P2 * (K2a + K2b)**3)

# Calculate the sum of velocity amplitudes for each system
K1_total = K1a + K1b
K2_total = K2a + K2b

# Calculate the mass ratio
# The units (years for P, km/s for K) are consistent between the two systems,
# so they will cancel out in the ratio, and no conversion to SI is needed.
try:
    mass_ratio = (P1 * (K1_total**3)) / (P2 * (K2_total**3))
except ZeroDivisionError:
    # This should not happen with the given data, but it's good practice
    print("Error: Division by zero. Check the input values for system 2.")
    mass_ratio = None

# The expected answer is D, which corresponds to ~0.4.
# The provided reasoning states the calculated value is 0.432.
expected_value_from_reasoning = 0.432
selected_option = 'D'
option_value = 0.4

# Check if the calculated ratio is close to the expected value and matches the selected option.
# We use a tolerance for floating-point comparison.
tolerance = 0.01

if mass_ratio is not None:
    # Check 1: Does our calculated value match the value from the reasoning?
    if not math.isclose(mass_ratio, expected_value_from_reasoning, rel_tol=tolerance):
        print(f"Incorrect. The reasoning states the calculated ratio is {expected_value_from_reasoning}, but the code calculates it as {mass_ratio:.3f}. There is a discrepancy in the calculation.")
    # Check 2: Does the calculated value correspond to the selected option?
    elif not math.isclose(mass_ratio, option_value, abs_tol=0.05): # Use a slightly larger tolerance for matching the option
        print(f"Incorrect. The calculated ratio is {mass_ratio:.3f}. While the reasoning might be correct, this value does not closely match option {selected_option} (~{option_value}).")
    else:
        # If both checks pass, the answer is correct.
        print("Correct")
