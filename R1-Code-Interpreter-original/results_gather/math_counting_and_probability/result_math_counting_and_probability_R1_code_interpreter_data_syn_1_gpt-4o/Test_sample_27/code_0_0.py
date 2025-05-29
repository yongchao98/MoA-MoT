# Given values
binom_19_12 = 50388
binom_17_11 = 12376

# Calculate binom(17, 5) using the relation with binom(19, 12)
binom_17_5 = binom_19_12 * 12 // (19 - 12 + 1)

# Calculate binom(18, 12) using Pascal's identity
binom_18_12 = binom_17_11 + binom_17_5

print(binom_18_12)