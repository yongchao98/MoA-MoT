# Statement III: The cumulative duration of the three shortest U.S. presidential
# administrations could fit within a local year at this location.

# The orbital period of 16 Cygni Bb in Earth days.
orbital_period_16_cygni_bb = 799.5

# The three shortest U.S. presidential administrations by days served.
# William Henry Harrison: March 4, 1841 – April 4, 1841
harrison_term = 31
# James A. Garfield: March 4, 1881 – September 19, 1881
garfield_term = 199
# Zachary Taylor: March 4, 1849 – July 9, 1850
taylor_term = 492

# Calculate the sum of these terms
sum_of_terms = harrison_term + garfield_term + taylor_term

# Print the values and the comparison
print("Analysis for Statement III:")
print(f"Orbital Period of 16 Cygni Bb: {orbital_period_16_cygni_bb} days")
print("\nThree Shortest U.S. Presidential Administrations:")
print(f"1. William Henry Harrison: {harrison_term} days")
print(f"2. James A. Garfield: {garfield_term} days")
print(f"3. Zachary Taylor: {taylor_term} days")

print(f"\nSum of shortest administrations: {harrison_term} + {garfield_term} + {taylor_term} = {sum_of_terms} days")

# Compare the sum with the orbital period
comparison = sum_of_terms < orbital_period_16_cygni_bb
print(f"\nCan the cumulative duration ({sum_of_terms} days) fit within a local year ({orbital_period_16_cygni_bb} days)?")
print(f"Result: {comparison} (True, the statement is correct)")
