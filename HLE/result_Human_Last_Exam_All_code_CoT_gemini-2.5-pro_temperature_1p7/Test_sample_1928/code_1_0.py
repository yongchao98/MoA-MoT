# Statement III: The cumulative duration of the three shortest U.S. presidential administrations
# could fit within a local year at this location.

# Orbital period of 16 Cygni Bb in days (a "local year").
orbital_period_16_cygni_bb = 799.5

# Durations of the three shortest U.S. presidential administrations in days.
# William Henry Harrison (March 4, 1841 – April 4, 1841)
harrison_term = 31
# James A. Garfield (March 4, 1881 – September 19, 1881)
garfield_term = 199
# Zachary Taylor (March 4, 1849 – July 9, 1850)
taylor_term = 492

# Calculate the cumulative duration.
cumulative_duration = harrison_term + garfield_term + taylor_term

# Print the calculation and the result.
print("Step 1: Calculate the sum of the three shortest presidential terms.")
print(f"Equation: {harrison_term} days (Harrison) + {garfield_term} days (Garfield) + {taylor_term} days (Taylor) = {cumulative_duration} days")
print("\nStep 2: Compare the sum to 16 Cygni Bb's orbital period.")
print(f"Cumulative duration of administrations: {cumulative_duration} days")
print(f"Orbital period of 16 Cygni Bb: {orbital_period_16_cygni_bb} days")

# Check if the cumulative duration fits within the orbital period.
is_statement_true = cumulative_duration < orbital_period_16_cygni_bb
print(f"\nIs the cumulative duration less than the orbital period? {is_statement_true}")

if is_statement_true:
    print("Conclusion: Statement III is true because the combined 722 days of the three shortest presidencies is less than the 799.5-day year on 16 Cygni Bb.")
else:
    print("Conclusion: Statement III is false.")
