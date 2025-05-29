# Barbara's circle
barbara_female_colleagues = 1  # Matilda is the female colleague of Barbara

# Charles's circle
charles_female_colleagues = 5  # Total female colleagues in Charles's circle
matilda_female_colleagues_in_charles_circle = charles_female_colleagues - 1  # Excluding Matilda herself

# Total female colleagues Matilda has
total_female_colleagues = barbara_female_colleagues + matilda_female_colleagues_in_charles_circle

print(total_female_colleagues)