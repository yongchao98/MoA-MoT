# Step 1 & 2: Define the historical figures.
# By 1823, the population of the Talbot Settlement was estimated at 12,000.
# The original land grant in 1803 was for 5,000 acres.
# Colonel Talbot eventually controlled the settlement of over 650,000 acres.
settlers_by_1823 = 12000
original_grant_acres = 5000
total_claimed_acres = 650000

# Step 3: Calculate the difference in acreage.
acreage_difference = total_claimed_acres - original_grant_acres

# Step 4: Print the results.
print(f"Number of destitute migrants settled by 1823: {settlers_by_1823}")
print("\nTo find how much larger the claimed acreage was than the original grant, we subtract the original grant from the total claimed acreage.")
print(f"The equation is: {total_claimed_acres} - {original_grant_acres} = {acreage_difference}")
print(f"\nThe acreage he claimed was {acreage_difference} acres larger than the original grant.")