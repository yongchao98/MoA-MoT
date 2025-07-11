# First, you may need to install the 'treys' library.
# You can do this by running the following command in your terminal or command prompt:
# pip install treys

try:
    from treys import Card, Evaluator, Deck
except ImportError:
    print("The 'treys' library is not installed.")
    print("Please install it by running: pip install treys")
    exit()

def simulate_equity(hero_hand_str, villain_hand_str, num_simulations=100000):
    """
    Simulates the equity of the hero's hand against the villain's hand.
    """
    hero_hand = [Card.new(c) for c in hero_hand_str]
    villain_hand = [Card.new(c) for c in villain_hand_str]
    evaluator = Evaluator()

    hero_wins = 0
    ties = 0

    # Create a deck once and store the integer representations of hole cards
    hero_hole_cards = [Card.new(hero_hand_str[0]), Card.new(hero_hand_str[1])]
    villain_hole_cards = [Card.new(villain_hand_str[0]), Card.new(villain_hand_str[1])]
    
    all_hole_cards = hero_hole_cards + villain_hole_cards

    for _ in range(num_simulations):
        deck = Deck()
        
        # Remove hole cards from the deck
        for card in all_hole_cards:
            deck.cards.remove(card)

        board = deck.draw(5)

        # The treys library uses integer representations for speed.
        # Lower score is better.
        hero_score = evaluator.evaluate(hero_hole_cards, board)
        villain_score = evaluator.evaluate(villain_hole_cards, board)

        if hero_score < villain_score:
            hero_wins += 1
        elif hero_score == villain_score:
            ties += 1

    equity = (hero_wins + ties / 2) / num_simulations
    return equity

# --- Define Hands ---
# Our hand: Two black aces
hero_hand = ['As', 'Ac']

# Opponent hands: Three different red suited hands (we'll use hearts)
villain_qj = ['Qh', 'Jh']
villain_qt = ['Qh', 'Th']
villain_q9 = ['Qh', '9h']

# --- Run Simulations ---
print("Running simulations to calculate equity for Aces...")
print("This may take a minute.\n")

equity_vs_qj = simulate_equity(hero_hand, villain_qj)
equity_vs_qt = simulate_equity(hero_hand, villain_qt)
equity_vs_q9 = simulate_equity(hero_hand, villain_q9)

# --- Display Results ---
print("--- Equity Calculation Results ---")
print(f"Equity for {hero_hand[0]}{hero_hand[1]} vs. {villain_qj[0]}{villain_qj[1]} (QJ suited): {equity_vs_qj:.2%}")
print(f"Equity for {hero_hand[0]}{hero_hand[1]} vs. {villain_qt[0]}{villain_qt[1]} (QT suited): {equity_vs_qt:.2%}")
print(f"Equity for {hero_hand[0]}{hero_hand[1]} vs. {villain_q9[0]}{villain_q9[1]} (Q9 suited): {equity_vs_q9:.2%}")
print("\n--- Conclusion ---")
print("We want to find the hand that minimizes our equity.")

equities = {
    "QJ suited": equity_vs_qj,
    "QT suited": equity_vs_qt,
    "Q9 suited": equity_vs_q9,
}

# Find the opponent hand that results in the lowest equity for us
worst_case_scenario = min(equities, key=equities.get)
min_equity_value = equities[worst_case_scenario]

print(f"\nThe lowest equity for our Aces ({min_equity_value:.2%}) occurs against {worst_case_scenario}.")
print("This is because suited connectors (like QJ) have a higher potential to make straights and flushes compared to gapped hands (like QT or Q9).")
print("Therefore, QJ suited is the hand you least want to see.")

<<<A>>>