# Define the rounds and results
rounds = [
    ('A', 'draw'),
    ('D', 'draw'),
    ('B', 'loss'),
    ('C', 'loss')
]

# Possible cards
cards = ['A', 'B', 'C', 'D']

# Initialize opponent's cards
opponent_cards = [None] * 4

# Deduce opponent's cards based on the rules
for i, (my_card, result) in enumerate(rounds):
    if result == 'draw':
        if my_card == 'A':
            opponent_cards[i] = 'D'  # Worm (D) is one step apart from Stick (A)
        elif my_card == 'D':
            opponent_cards[i] = 'C'  # Chicken (C) is one step apart from Worm (D)
    elif result == 'loss':
        if my_card == 'B':
            opponent_cards[i] = 'A'  # Stick (A) beats Tiger (B)
        elif my_card == 'C':
            opponent_cards[i] = 'B'  # Tiger (B) beats Chicken (C)

print(opponent_cards)