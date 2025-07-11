# Historical data for the calculation
settlers_by_1823 = 12000
original_grant_acres = 5000
final_claimed_acres = 650000

# Calculate the difference in acreage
acreage_difference = final_claimed_acres - original_grant_acres

# Print the answers to both parts of the question
print(f"Between 1803 and 1823, a total of {settlers_by_1823} migrants settled in the Talbot Settlement.")
print(f"The acreage Colonel Talbot ultimately claimed was {acreage_difference:,} acres larger than his original land grant.")
print(f"This is calculated as: {final_claimed_acres:,} acres - {original_grant_acres:,} acres = {acreage_difference:,} acres.")
