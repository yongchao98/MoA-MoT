# Initial cards
cards = list("QOBSYIAVCJNYK")

# Define the rounds with the card to add and the card to remove
rounds = [
    ('P', 'C'),
    ('U', 'J'),
    ('W', 'S'),
    ('Z', 'Y'),
    ('[', None)  # In the final round, we don't remove a card
]

# Function to check for "Peng"
def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

# Function to check for "Chi"
def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Process each round
for new_card, remove_card in rounds:
    # Check for "Peng"
    if check_peng(cards, new_card):
        result = 1
    # Check for "Chi"
    elif check_chi(cards, new_card):
        result = 2
    else:
        result = 0
    
    # Remove the specified card if it's not the final round
    if remove_card:
        cards.remove(remove_card)
    
    # Add the new card
    cards.append(new_card)

# Print the result of the final round
print(result)