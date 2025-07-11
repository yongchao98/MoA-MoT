# Data based on historical records for Colonel Thomas Talbot's settlement.
settlers_by_1823 = 20000
original_land_grant_acres = 5000
total_claimed_acres = 650000

# Calculate the difference in acreage
acreage_difference = total_claimed_acres - original_land_grant_acres

# Print the answers to the user's questions
print(f"Number of destitute migrants settled between 1803 and 1823: {settlers_by_1823}")
print(f"The acreage Colonel Talbot claimed was larger than the original land grant by {total_claimed_acres} - {original_land_grant_acres} acres.")
print(f"This is a difference of {acreage_difference} acres.")