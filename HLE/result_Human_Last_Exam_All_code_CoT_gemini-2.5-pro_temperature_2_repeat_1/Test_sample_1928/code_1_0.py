import math

# --- Data for Calculation ---

# Orbital period of 16 Cygni Bb in days
orbital_period_16_cyg_bb = 799.5

# Durations of the three shortest U.S. presidential administrations in days
# William Henry Harrison (March 4, 1841 – April 4, 1841)
harrison_term = 31
# James A. Garfield (March 4, 1881 – September 19, 1881)
garfield_term = 199
# Zachary Taylor (March 4, 1849 – July 9, 1850)
# 1 year (365 days) + 127 days in 1850
taylor_term = 492

# --- Calculation ---

# Sum the durations of the administrations
cumulative_duration = harrison_term + garfield_term + taylor_term

# --- Output ---

print("This script verifies statement III.")
print("-" * 35)
print(f"Orbital period of 16 Cygni Bb: {orbital_period_16_cyg_bb} days")
print("\nShortest U.S. Presidential Administrations:")
print(f"1. William Henry Harrison: {harrison_term} days")
print(f"2. James A. Garfield:     {garfield_term} days")
print(f"3. Zachary Taylor:        {taylor_term} days")
print("\nCumulative duration calculation:")

# Show the full equation as requested
print(f"{harrison_term} + {garfield_term} + {taylor_term} = {cumulative_duration} days")

# Final comparison
is_shorter = cumulative_duration < orbital_period_16_cyg_bb
print(f"\nIs the cumulative duration ({cumulative_duration} days) less than a year on 16 Cygni Bb ({orbital_period_16_cyg_bb} days)?")
print(f"Result: {is_shorter}")

if is_shorter:
    print("\nConclusion: Statement III is true.")
else:
    print("\nConclusion: Statement III is false.")
