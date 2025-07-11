import math

# Part 1: Data regarding settlers
start_year = 1803
end_year = 1823
settler_population_by_1823 = 20000

# Part 2: Data regarding acreage
original_grant_acres = 5000
final_claimed_acres = 650000

# Calculate the difference in acreage
acreage_difference = final_claimed_acres - original_grant_acres

# Print the answers
print(f"Between {start_year} and {end_year}, the Talbot Settlement grew to a population of approximately {settler_population_by_1823} migrants.")
print("\n---")
print(f"The acreage Colonel Talbot ultimately claimed was {acreage_difference:,} acres larger than his original land grant.")
print(f"Calculation: {final_claimed_acres:,} (claimed) - {original_grant_acres:,} (original grant) = {acreage_difference:,} acres")