# Number of petals per orchid and daisy
petals_per_orchid = 5
petals_per_daisy = 10

# Number of orchids and daisies
num_orchids = 5
num_daisies = 8516550

# Total petals for orchids and daisies
total_petals_orchids = num_orchids * petals_per_orchid
total_petals_daisies = num_daisies * petals_per_daisy

# Difference in petals
difference = total_petals_daisies - total_petals_orchids

print(difference)