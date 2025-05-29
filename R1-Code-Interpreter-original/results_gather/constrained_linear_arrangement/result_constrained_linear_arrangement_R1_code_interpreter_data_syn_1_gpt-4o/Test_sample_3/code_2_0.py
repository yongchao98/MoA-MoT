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
    all_cards = {'A', 'B', 'C', 'D'}
    
    for i, (your_card, result) in enumerate(zip(your_cards, results)):
        if result == 'draw':
            opponent_card = your_card
        elif result == 'loss':
            opponent_card = rules[your_card]
        
        # Ensure the card hasn't been used already
        if opponent_card in used_cards:
            # Find an unused card that results in a draw or loss
            for card in all_cards - used_cards:
                if (result == 'draw' and card == your_card) or (result == 'loss' and rules[your_card] == card):
                    opponent_card = card
                    break
        
        opponent_cards.append(opponent_card)
        used_cards.add(opponent_card)
    
    return opponent_cards

# Determine the opponent's cards
opponent_cards = determine_opponent_cards(your_cards, results)
print(opponent_cards)