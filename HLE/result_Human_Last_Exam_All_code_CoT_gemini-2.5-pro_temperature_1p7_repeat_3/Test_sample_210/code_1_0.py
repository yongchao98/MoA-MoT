# Part 1: Find the number of settlers
# Historical sources indicate that by the mid-1820s, the population
# of the Talbot Settlement was approximately 20,000.
settlers_1823 = 20000

# Part 2: Calculate the difference in acreage
# The original land grant was 5,000 acres.
original_grant_acres = 5000
# The total acreage he eventually controlled was about 650,000 acres.
total_claimed_acres = 650000

# Calculate the difference
acreage_difference = total_claimed_acres - original_grant_acres

# Print the results
print(f"Number of destitute migrants settled between 1803 and 1823: {settlers_1823}")
print("\nCalculation for the difference in acreage:")
print(f"The acreage Colonel Talbot claimed ({total_claimed_acres} acres) was larger than the original grant ({original_grant_acres} acres) by:")
print(f"{total_claimed_acres} - {original_grant_acres} = {acreage_difference} acres")