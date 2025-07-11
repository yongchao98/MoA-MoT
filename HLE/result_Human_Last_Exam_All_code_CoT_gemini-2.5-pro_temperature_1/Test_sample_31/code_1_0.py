import math

# Step 1: Calculate the theoretical minimum number of stickers needed.

num_pll_cases = 21
num_sticker_colors = 4

print("To solve this, we first establish a theoretical minimum number of stickers.")
print("We need to distinguish between {} different PLL cases.".format(num_pll_cases))
print("Each non-top-facing sticker can be one of {} possible colors.".format(num_sticker_colors))
print("Let 'n' be the number of stickers we observe. The number of possible outcomes is 4^n.")
print("We need to find the smallest 'n' where 4^n is greater than or equal to 21.\n")

# Check n=2
n2 = 2
outcomes_n2 = num_sticker_colors**n2
print("For n = {}, the equation is {} ^ {} = {}.".format(n2, num_sticker_colors, n2, outcomes_n2))
print("Since {} is less than {}, 2 stickers are not enough.\n".format(outcomes_n2, num_pll_cases))

# Check n=3
n3 = 3
outcomes_n3 = num_sticker_colors**n3
print("For n = {}, the equation is {} ^ {} = {}.".format(n3, num_sticker_colors, n3, outcomes_n3))
print("Since {} is greater than {}, the theoretical minimum number of stickers is 3.\n".format(outcomes_n3, num_pll_cases))

# Step 2 & 3: Analyze the sufficiency of the theoretical minimum.
print("However, this theoretical minimum is not practically sufficient.")
print("No matter which 3 sticker positions you choose, there will always be at least two PLL cases that look identical.")
print("For example, if you observe 3 stickers on pieces that are not affected by a U-perm (a 3-cycle of edges), you cannot distinguish that U-perm from a solved state.")
print("This means 3 stickers are not enough to resolve all ambiguities.\n")
print("Therefore, we must look at one more sticker.")
print("The minimum number of stickers that must be seen to fully identify any PLL case is 4.")
