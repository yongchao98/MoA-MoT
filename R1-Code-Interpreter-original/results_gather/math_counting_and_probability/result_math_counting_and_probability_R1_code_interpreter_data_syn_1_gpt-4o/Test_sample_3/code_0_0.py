# Given values
binom_15_6 = 5005
binom_15_8 = 6435

# Using symmetry property
binom_15_7 = binom_15_8
binom_15_9 = binom_15_6

# Calculate binom_16_8 and binom_16_9
binom_16_8 = binom_15_7 + binom_15_8
binom_16_9 = binom_15_8 + binom_15_9

# Calculate binom_17_9
binom_17_9 = binom_16_8 + binom_16_9

print(binom_17_9)