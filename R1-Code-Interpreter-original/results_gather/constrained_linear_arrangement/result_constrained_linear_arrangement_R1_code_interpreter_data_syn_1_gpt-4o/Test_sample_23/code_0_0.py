# Deduce opponent's cards based on the given rounds and results
rounds = [
    ('A', 'loss'),  # Round 1
    ('B', 'draw'),  # Round 2
    ('D', 'win'),   # Round 3
    ('C', 'draw')   # Round 4
]

# Possible cards
cards = ['A', 'B', 'C', 'D']

# Initialize opponent's cards
opponent_cards = [None] * 4

# Deduce each round
for i, (my_card, result) in enumerate(rounds):
    if result == 'loss':
        if my_card == 'A':
            opponent_cards[i] = 'C'  # Stick loses to Chicken
        elif my_card == 'B':
            opponent_cards[i] = 'D'  # Tiger loses to Worm
        elif my_card == 'C':
            opponent_cards[i] = 'A'  # Chicken loses to Stick
        elif my_card == 'D':
            opponent_cards[i] = 'B'  # Worm loses to Tiger
    elif result == 'draw':
        opponent_cards[i] = my_card  # Draw means same card
    elif result == 'win':
        if my_card == 'A':
            opponent_cards[i] = 'B'  # Stick beats Tiger
        elif my_card == 'B':
            opponent_cards[i] = 'C'  # Tiger beats Chicken
        elif my_card == 'C':
            opponent_cards[i] = 'D'  # Chicken beats Worm
        elif my_card == 'D':
            opponent_cards[i] = 'A'  # Worm beats Stick

# Print the deduced opponent's cards
print(opponent_cards)