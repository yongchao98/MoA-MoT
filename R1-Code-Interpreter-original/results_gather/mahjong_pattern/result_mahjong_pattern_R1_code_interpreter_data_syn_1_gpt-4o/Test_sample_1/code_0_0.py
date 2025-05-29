# Initial cards
initial_cards = "JHUFWIXNAPYRX"

# New card added
new_card = "A"

# Count the occurrences of the new card in the initial set
count_new_card = initial_cards.count(new_card)

# Determine the result based on the count
result = 1 if count_new_card >= 2 else 0

# Output the result
print(result)