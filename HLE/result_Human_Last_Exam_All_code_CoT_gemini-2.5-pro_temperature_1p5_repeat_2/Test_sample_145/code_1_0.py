import sys
import random

# This script requires the 'deuces' library.
# If you don't have it, please install it by running: pip install deuces
try:
    from deuces import Card, Evaluator, Deck
except ImportError:
    print("Error: The 'deuces' library is not installed.")
    print("Please install it to run this script by executing: pip install deuces")
    sys.exit()

def calculate_equity(hero_cards_str, villain_cards_str, num_simulations=100000):
    """
    Calculates the equity of the hero's hand against the villain's hand
    through Monte Carlo simulation.
    """
    evaluator = Evaluator()
    hero_hand = [Card.new(c) for c in hero_cards_str]
    villain_hand = [Card.new(c) for c in villain_cards_str]
    
    hero_wins = 0
    ties = 0

    # Run the simulation
    for _ in range(num_simulations):
        deck = Deck()
        
        # Remove known cards from the deck
        deck.cards.remove(hero_hand[0])
        deck.cards.remove(hero_hand[1])
        deck.cards.remove(villain_hand[0])
        deck.cards.remove(villain_hand[1])
        
        board = deck.draw(5)
        
        # Evaluate hand strengths (lower score is better in 'deuces')
        hero_score = evaluator.evaluate(hero_hand, board)
        villain_score = evaluator.evaluate(villain_hand, board)

        if hero_score < villain_score:
            hero_wins += 1
        elif hero_score == villain_score:
            ties += 1

    # Calculate equity percentage
    equity = (hero_wins + ties / 2.0) / num_simulations * 100
    return equity

# --- Main Logic ---

# Define the hands
# We have two black aces, e.g., Ace of Spades and Ace of Clubs
hero_hand = ['As', 'Ac']

# The opponent has one of three red suited hands
villain_hands = {
    "A. QJ": ['Qh', 'Jh'],
    "B. QT": ['Qh', 'Th'],
    "C. Q9": ['Qh', '9h'],
}

print("Running poker equity simulations...")
print(f"Our Hand: {Card.get_pretty_str(Card.new(hero_hand[0]))}, {Card.get_pretty_str(Card.new(hero_hand[1]))}")
print("-" * 45)

results = {}
# Run simulation for each opponent hand
for name, hand in villain_hands.items():
    equity = calculate_equity(hero_hand, hand)
    results[name] = equity
    pretty_hand = f"{Card.get_pretty_str(Card.new(hand[0]))}, {Card.get_pretty_str(Card.new(hand[1]))}"
    print(f"Our equity vs. Opponent's {name} ({pretty_hand}): {equity:.2f}%")

# Determine which hand we least like to see (the one that minimizes our equity)
min_equity_hand = min(results, key=results.get)

print("-" * 45)
print(f"The opponent hand that results in the lowest equity for our Aces is {min_equity_hand}.")
print("Therefore, this is the hand we least like to see.")

# Map the result to the provided answer choices
final_answer_letter = min_equity_hand.split('.')[0]
print(f"\nFinal Answer based on the simulation is {min_equity_hand}.")

<<<A>>>