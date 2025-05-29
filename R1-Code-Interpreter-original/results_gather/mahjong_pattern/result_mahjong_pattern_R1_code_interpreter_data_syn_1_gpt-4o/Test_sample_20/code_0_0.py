# Initial set of cards
cards = list("VUINWATZHPRPQ")

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

# Rounds information
rounds = [
    ('P', 'A'),
    ('J', 'R'),
    ('W', 'J'),
    ('P', 'T'),
    ('J', None)  # No card to remove in the final round
]

# Process each round
for new_card, remove_card in rounds:
    if check_peng(cards, new_card):
        result = 1  # Peng
    elif check_chi(cards, new_card):
        result = 2  # Chi
    else:
        result = 0  # Pass
    
    # Add the new card
    cards.append(new_card)
    
    # Remove the specified card if any
    if remove_card:
        cards.remove(remove_card)

# Output the result of the final round
print(result)