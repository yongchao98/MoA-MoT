# Define the cards and results
rounds = [
    ('C', 'draw'),
    ('A', 'win'),
    ('D', 'loss'),
    ('B', 'draw')
]

# Define the rules
rules = {
    'A': 'B',  # Stick beats Tiger
    'B': 'C',  # Tiger beats Chicken
    'C': 'D',  # Chicken beats Worm
    'D': 'A'   # Worm beats Stick
}

# Deduce opponent's cards
opponent_cards = []

for i, (my_card, result) in enumerate(rounds):
    if result == 'draw':
        if my_card == 'A':
            opponent_cards.append('D')
        elif my_card == 'B':
            opponent_cards.append('A')
        elif my_card == 'C':
            opponent_cards.append('C')
        elif my_card == 'D':
            opponent_cards.append('B')
    elif result == 'win':
        opponent_cards.append(rules[my_card])
    elif result == 'loss':
        for card, beats in rules.items():
            if beats == my_card:
                opponent_cards.append(card)
                break

print(opponent_cards)