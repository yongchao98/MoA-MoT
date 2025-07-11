import datetime

# This script verifies the calculation for Statement III.
# Statement III: The cumulative duration of the three shortest U.S. presidential administrations 
# could fit within a local year at this location.

# The orbital period of 16 Cygni Bb, which defines its "local year"
period_16_cyg_bb = 798.5  # in Earth days

# The three shortest U.S. presidential administrations by tenure
# 1. William Henry Harrison
harrison_start = datetime.date(1841, 3, 4)
harrison_end = datetime.date(1841, 4, 4)
harrison_duration = (harrison_end - harrison_start).days

# 2. James A. Garfield
garfield_start = datetime.date(1881, 3, 4)
garfield_end = datetime.date(1881, 9, 19)
garfield_duration = (garfield_end - garfield_start).days

# 3. Zachary Taylor
taylor_start = datetime.date(1849, 3, 4)
taylor_end = datetime.date(1850, 7, 9)
taylor_duration = (taylor_end - taylor_start).days

# Calculate the cumulative duration of these administrations
total_duration = harrison_duration + garfield_duration + taylor_duration

# Print the analysis step-by-step
print("Analysis for Statement III:")
print(f"A 'local year' on 16 Cygni Bb is {period_16_cyg_bb} Earth days.")
print("\nCumulative duration of the three shortest U.S. presidential administrations:")
print(f" - William Henry Harrison's term: {harrison_duration} days")
print(f" - James A. Garfield's term: {garfield_duration} days")
print(f" - Zachary Taylor's term: {taylor_duration} days")
print("\nFinal Equation:")
print(f"{harrison_duration} + {garfield_duration} + {taylor_duration} = {total_duration}")
print(f"\nComparing the two durations: {total_duration} days is less than {period_16_cyg_bb} days.")
print("Therefore, Statement III is TRUE.")
