# Initial cards
cards = list("PCHWEHUWNAGXO")

# Define the rounds
rounds = [
    ('P', 'W'),  # Round 1
    ('Z', 'N'),  # Round 2
    ('I', 'H'),  # Round 3
    ('H', 'H'),  # Round 4
    ('P', None)  # Round 5, no removal needed
]

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

# Process each round
for new_card, remove_card in rounds:
    cards.append(new_card)
    if remove_card:
        cards.remove(remove_card)
    
    # Check for "Peng"
    if check_peng(cards, new_card):
        result = 1
    # Check for "Chi"
    elif check_chi(cards, new_card):
        result = 2
    else:
        result = 0

# Print the result of the final round
print(result)