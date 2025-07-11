# Historical data for the Talbot Settlement
settlers_by_1823 = 20000
original_grant_acres = 5000
total_claimed_acres = 650000

# Calculate the difference in acreage
acreage_difference = total_claimed_acres - original_grant_acres

# Print the answers
print(f"Number of destitute migrants settled between 1803 and 1823: {settlers_by_1823}")
print("\nCalculation for the difference in acreage:")
print(f"Total Acreage Claimed: {total_claimed_acres}")
print(f"Original Land Grant: {original_grant_acres}")
print(f"The acreage Colonel Talbot claimed was larger than the original grant by: {total_claimed_acres} - {original_grant_acres} = {acreage_difference} acres.")

# Final answer in the requested format
final_answer = (
    f"Between 1803 and 1823, {settlers_by_1823} destitute migrants settled as a result of the land grant. "
    f"The acreage Colonel Talbot eventually claimed was {acreage_difference} acres larger than the original "
    f"{original_grant_acres} acre grant."
)

# This final line is for the format requirement and will not be visible in a normal run.
# print(f"\n<<<{final_answer}>>>")