# Sarah's circle
sarah_male_colleagues = 2
sarah_female_colleagues = 4

# Thomas's circle
thomas_male_colleagues = 1
thomas_female_colleagues = 2

# Matilda is common in both circles
# Total unique female colleagues for Matilda
unique_female_colleagues = sarah_female_colleagues + thomas_female_colleagues - 1  # Subtract 1 for Matilda herself

print(unique_female_colleagues)