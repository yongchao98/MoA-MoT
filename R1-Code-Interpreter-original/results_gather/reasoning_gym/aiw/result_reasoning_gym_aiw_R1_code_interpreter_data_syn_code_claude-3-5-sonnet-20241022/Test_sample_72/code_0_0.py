# Linda's circle (including Linda)
linda_males = 4
linda_females = 3  # excluding Linda herself
linda_total = linda_males + linda_females + 1  # +1 for Linda

# Michael's circle (including Michael)
michael_males = 1  # excluding Michael himself
michael_females = 6
michael_total = michael_males + michael_females + 1  # +1 for Michael

# We know Matilda is in both circles
# Let's calculate Matilda's female colleagues:

# In Linda's circle:
# - Linda herself
# - The other 2 females (3 total - Matilda)
linda_circle_females_for_matilda = 1 + (linda_females - 1)

# In Michael's circle:
# - The other 5 females (6 total - Matilda)
michael_circle_females_for_matilda = michael_females - 1

# Total unique female colleagues of Matilda
total_female_colleagues = linda_circle_females_for_matilda + michael_circle_females_for_matilda

print(f"Matilda's female colleagues: {total_female_colleagues}")