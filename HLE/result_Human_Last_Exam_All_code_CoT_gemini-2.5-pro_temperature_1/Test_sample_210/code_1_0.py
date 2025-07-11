# Historical data for Colonel Thomas Talbot's Settlement
settlers_by_1823 = 20000
original_land_grant_acres = 5000
final_claimed_acres = 65000

# Calculate the difference in acreage
acreage_difference = final_claimed_acres - original_land_grant_acres

# Print the answers to the user's questions
print(f"Number of migrants settled between 1803 and 1823: {settlers_by_1823}")
print("\nCalculation for the difference in acreage:")
print(f"The acreage Talbot eventually claimed was {final_claimed_acres - original_land_grant_acres} acres larger than the original grant.")
print(f"Equation: {final_claimed_acres} (claimed) - {original_land_grant_acres} (original) = {acreage_difference} acres")
