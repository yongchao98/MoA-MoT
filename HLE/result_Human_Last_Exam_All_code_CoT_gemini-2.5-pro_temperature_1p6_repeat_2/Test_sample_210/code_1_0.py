# Step 1: Define variables based on historical data.
# By the early 1820s, the population of the Talbot Settlement had grown significantly.
migrants_settled_by_1823 = 20000

# Colonel Talbot's original land grant in 1803.
original_land_grant_acres = 5000

# The total acreage Talbot eventually claimed for his personal estate.
final_claimed_acres = 65000

# Step 2: Calculate how much larger the final claimed acreage was.
acreage_difference = final_claimed_acres - original_land_grant_acres

# Step 3: Print the results to answer the user's questions.
print(f"As a result of the land grant, approximately {migrants_settled_by_1823} migrants had settled in the Talbot Settlement by 1823.")
print("\nTo determine how much larger the acreage he claimed was compared to the original grant, we perform the following calculation:")
print(f"Final Claimed Acres - Original Grant Acres = Increase in Acres")
print(f"{final_claimed_acres} - {original_land_grant_acres} = {acreage_difference}")
print(f"\nColonel Talbot's claimed acreage was {acreage_difference} acres larger than his original land grant.")
