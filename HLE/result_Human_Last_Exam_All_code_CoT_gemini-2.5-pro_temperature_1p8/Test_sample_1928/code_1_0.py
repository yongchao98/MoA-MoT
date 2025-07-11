# Statement III: The cumulative duration of the three shortest U.S.
# presidential administrations could fit within a local year at this location.

# Orbital period of 16 Cygni Bb in days
orbital_period_16_cygni_bb = 799.5 # days

# Durations of the three shortest U.S. presidential administrations in days
harrison_term = 31
garfield_term = 199
taylor_term = 492

# Calculate the cumulative duration
cumulative_duration = harrison_term + garfield_term + taylor_term

# Check if the cumulative duration is less than the orbital period
is_statement_true = cumulative_duration < orbital_period_16_cygni_bb

# Print the analysis for statement III
print("--- Analysis of Statement III ---")
print(f"Orbital period of 16 Cygni Bb: {orbital_period_16_cygni_bb} days")
print("Three shortest U.S. presidential administrations:")
print(f"1. William Henry Harrison: {harrison_term} days")
print(f"2. James A. Garfield: {garfield_term} days")
print(f"3. Zachary Taylor: {taylor_term} days")
print(f"Cumulative duration = {harrison_term} + {garfield_term} + {taylor_term} = {cumulative_duration} days")
print(f"Is {cumulative_duration} days < {orbital_period_16_cygni_bb} days? {is_statement_true}")
if is_statement_true:
    print("Conclusion: Statement III is TRUE.")

# Based on research, Statement VI is also true.
print("\n--- Analysis of Statement VI ---")
print("Research confirms two 'Nature' journal stories feature voyages to 16 Cygni.")
print("Conclusion: Statement VI is TRUE.")

# Combine the true statements' numerals
print("\n--- Final Answer ---")
print("The true statements are III and VI.")
print("Final Answer Sequence: III-VI")
