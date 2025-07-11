# This script verifies the calculation for Statement III.

# The orbital period of 16 Cygni Bb defines its "local year".
orbital_period_16_cygni_bb = 799.5  # in days

# The durations of the three shortest U.S. presidential administrations.
harrison_term = 31  # William Henry Harrison
garfield_term = 199 # James A. Garfield
taylor_term = 492   # Zachary Taylor

# Calculate the sum of the three shortest administrations.
cumulative_duration = harrison_term + garfield_term + taylor_term

# Display the information and the calculation.
print("--- Verifying Statement III ---")
print(f"A 'local year' on 16 Cygni Bb (its orbital period) is {orbital_period_16_cygni_bb} days.")
print("\nThe cumulative duration of the three shortest U.S. presidential administrations is calculated as follows:")
print(f"William Henry Harrison's term: {harrison_term} days")
print(f"James A. Garfield's term: {garfield_term} days")
print(f"Zachary Taylor's term: {taylor_term} days")
print("\nEquation:")
print(f"{harrison_term} + {garfield_term} + {taylor_term} = {cumulative_duration}")
print(f"\nTotal cumulative duration: {cumulative_duration} days.")

# Compare the two durations.
is_shorter = cumulative_duration < orbital_period_16_cygni_bb

print(f"\nIs the cumulative duration ({cumulative_duration} days) less than the local year ({orbital_period_16_cygni_bb} days)? {is_shorter}")
if is_shorter:
    print("Conclusion: Statement III is true.")
else:
    print("Conclusion: Statement III is false.")

# Based on the full analysis, the true statements are III and VI.
print("\nFinal Answer Sequence: III-VI")