# Initial cards
cards = list("VFGJVKQKLDPTH")

# Function to check for "Peng"
def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

# Function to check for "Chi"
def check_chi(cards, new_card):
    # Sort the cards and check for consecutive sequence
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Simulate each round
rounds = [
    ('V', 'V'),
    ('K', 'K'),
    ('L', 'V'),
    ('L', 'H'),
    ('R', None)  # No card to remove in the final round
]

result = 0
for new_card, remove_card in rounds:
    if check_peng(cards, new_card):
        result = 1
    elif check_chi(cards, new_card):
        result = 2
    else:
        result = 0
    
    # Add the new card
    cards.append(new_card)
    
    # Remove the specified card
    if remove_card:
        cards.remove(remove_card)

# Output the result of the final round
print(result)