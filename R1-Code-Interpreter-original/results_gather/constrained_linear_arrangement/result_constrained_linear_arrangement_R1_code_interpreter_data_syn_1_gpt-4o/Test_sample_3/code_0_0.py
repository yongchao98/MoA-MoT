# Define the cards and results
your_cards = ['A', 'B', 'C', 'D']
results = ['draw', 'draw', 'loss', 'loss']

# Define the rules
rules = {
    'A': 'B',  # Stick beats Tiger
    'B': 'C',  # Tiger beats Chicken
    'C': 'D',  # Chicken beats Worm
    'D': 'A'   # Worm beats Stick
}

# Function to determine opponent's card
def determine_opponent_cards(your_cards, results):
    opponent_cards = []
    used_cards = set()
    
    for i, (your_card, result) in enumerate(zip(your_cards, results)):
        if result == 'draw':
            if your_card == 'A':
                opponent_card = 'A'
            elif your_card == 'B':
                opponent_card = 'B'
        elif result == 'loss':
            opponent_card = rules[your_card]
        
        opponent_cards.append(opponent_card)
        used_cards.add(opponent_card)
    
    return opponent_cards

# Determine the opponent's cards
opponent_cards = determine_opponent_cards(your_cards, results)
print(opponent_cards)