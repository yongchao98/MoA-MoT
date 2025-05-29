# Define the sets
cats = set(range(1, 11))  # Let's assume there are 10 cats
adults = set(range(1, 21))  # Let's assume there are 20 adults
humans = set(range(11, 21))  # Let's assume humans are adults numbered 11 to 20

# All cats are adults
assert cats.issubset(adults)

# Some adults are not humans
non_human_adults = adults - humans
assert len(non_human_adults) > 0

# Check if some cats are not humans
cats_not_humans = cats - humans
result = len(cats_not_humans) > 0

print(result)