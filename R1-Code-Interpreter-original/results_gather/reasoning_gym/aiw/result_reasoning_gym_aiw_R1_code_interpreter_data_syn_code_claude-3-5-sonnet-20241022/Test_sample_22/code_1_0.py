# Jennifer's circle females (excluding Matilda and Jennifer)
jennifer_circle_females = 3 - 1  # -1 for Matilda

# Joseph's circle females (excluding Matilda and Joseph)
joseph_circle_females = 4 - 1  # -1 for Matilda

# Total female colleagues of Matilda
# We add both numbers because these are distinct groups
total_female_colleagues = jennifer_circle_females + joseph_circle_females

print(f"Matilda has {total_female_colleagues} female colleagues")