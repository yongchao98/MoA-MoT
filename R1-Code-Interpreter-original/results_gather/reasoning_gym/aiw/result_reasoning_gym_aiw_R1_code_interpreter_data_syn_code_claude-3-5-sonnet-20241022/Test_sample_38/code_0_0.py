# Create sets for Elizabeth's and Joseph's circles
elizabeth_total = 6  # total colleagues
elizabeth_males = 4
elizabeth_females = 2

joseph_total = 4  # total colleagues
joseph_males = 3
joseph_females = 1

# We know Matilda is in both circles
# Let's find Matilda's female colleagues

# In Elizabeth's circle:
# Matilda is one of the 2 female colleagues, so there is 1 other female
elizabeth_female_colleagues_for_matilda = elizabeth_females - 1  # excluding Matilda herself

# In Joseph's circle:
# Matilda is the only female
joseph_female_colleagues_for_matilda = 0  # no other females in Joseph's circle

# Total female colleagues for Matilda
total_female_colleagues = elizabeth_female_colleagues_for_matilda + joseph_female_colleagues_for_matilda

print(f"Matilda has {total_female_colleagues} female colleague(s)")