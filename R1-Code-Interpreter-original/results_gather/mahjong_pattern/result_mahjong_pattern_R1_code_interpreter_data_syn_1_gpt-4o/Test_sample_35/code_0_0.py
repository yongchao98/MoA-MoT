# Initial cards
cards = list("GSNRTAVCJKZYF")

# Define the rounds
rounds = [
    ('J', 'K'),
    ('J', 'V'),
    ('N', 'A'),
    ('H', 'Z'),
    ('J', None)  # No card to remove in the final round
]

# Function to check for "Peng"
def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

# Function to check for "Chi"
def check_chi(cards, new_card):
    cards.append(new_card)
    cards.sort()
    for i in range(len(cards) - 2):
        if ord(cards[i+1]) == ord(cards[i]) + 1 and ord(cards[i+2]) == ord(cards[i]) + 2:
            return True
    return False

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
    
    # Remove the specified card
    if remove_card:
        cards.remove(remove_card)

# Output the result of the final round
print(result)