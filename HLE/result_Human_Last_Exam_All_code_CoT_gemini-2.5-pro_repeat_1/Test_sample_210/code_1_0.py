# Historical data for the Talbot Settlement
original_land_grant_acres = 5000
total_acreage_claimed = 650000
settlers_by_1823 = 12000

# Calculate the difference in acreage
acreage_difference = total_acreage_claimed - original_land_grant_acres

# Print the final answer, showing the calculation for the acreage
print(f"By 1823, {settlers_by_1823} settlers had migrated to the Talbot Settlement.")
print(f"The acreage he eventually claimed was {total_acreage_claimed:,} acres, which was {acreage_difference:,} acres larger than the original {original_land_grant_acres:,} acre grant ({total_acreage_claimed:,} - {original_land_grant_acres:,} = {acreage_difference:,}).")