# Patricia's circle (excluding Matilda)
patricia_males = 6
patricia_females = 5  # 6 minus Matilda

# William's circle (excluding Matilda)
william_males = 3
william_females = 5  # 6 minus Matilda

# Matilda's female colleagues would be the sum of females in both circles
# Since Matilda is in both circles, we need to exclude her from both counts
total_female_colleagues = patricia_females + william_females

# Since these are the only colleagues (no other connections outside these circles)
# and Matilda is the only overlap between circles
# The total number of female colleagues for Matilda is simply the sum

print(f"Matilda's female colleagues: {total_female_colleagues}")