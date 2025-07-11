import itertools

# First, I will explain the plan.
# The goal is to find out which of the three hands (QJ, QT, Q9 suited) is the strongest
# against a pair of black aces (AsAc). The strongest opponent hand is the one that gives
# the aces the lowest probability of winning (lowest equity).
#
# To do this, we can simulate the game by iterating through every possible 5-card board.
# 1. We start with a 52-card deck.
# 2. We remove our two aces (e.g., As, Ac) and the opponent's two cards (e.g., Qh, Jh).
# 3. This leaves 48 cards in the deck. We will check every combination of 5 cards from these 48.
# 4. For each combination (board), we evaluate our hand and the opponent's hand to see who wins.
# 5. We tally the wins, losses, and ties for our aces.
# 6. Equity is calculated as: (number of wins + number of ties / 2) / total number of boards.
#
# The script below uses the 'deuces' library to perform the hand evaluations.
# You may need to install it first by running: pip install deuces

try:
    from deuces import Card, Evaluator, Deck
except ImportError:
    print("This script requires the 'deuces' library.")
    print("Please run 'pip install deuces' in your terminal and then run this script again.")
    exit()

def calculate_equity_and_print(hero_cards_str, villain_cards_str):
    """
    Calculates the equity for a hero vs. villain matchup by enumerating all possible boards.
    It prints the detailed components of the equity calculation.
    """
    evaluator = Evaluator()
    hero_hand = [Card.new(c) for c in hero_cards_str]
    villain_hand = [Card.new(c) for c in villain_cards_str]

    # Create a deck and remove the known cards to get the 48 remaining cards
    deck = Deck()
    known_cards_str = set(hero_cards_str + villain_cards_str)
    deck_list = [c for c in deck.cards if Card.get_pretty_str(c) not in known_cards_str]
    
    hero_wins, ties, total_boards = 0, 0, 0

    # Iterate through all possible 5-card boards (C(48, 5) = 1,712,304 combinations)
    for board_tuple in itertools.combinations(deck_list, 5):
        total_boards += 1
        # The 'deuces' evaluator returns a numerical rank; a lower number is a better hand
        hero_score = evaluator.evaluate(hero_hand, list(board_tuple))
        villain_score = evaluator.evaluate(villain_hand, list(board_tuple))

        if hero_score < villain_score:
            hero_wins += 1
        elif hero_score == villain_score:
            ties += 1

    equity = (hero_wins + ties / 2) / total_boards
    
    hero_hand_pretty = f"{hero_cards_str[0]}{hero_cards_str[1]}"
    villain_hand_pretty = f"{villain_cards_str[0]}{villain_cards_str[1]}"
    
    print(f"Analysis for {hero_hand_pretty} vs. {villain_hand_pretty}:")
    print(f"Total unique boards: {total_boards}")
    print(f"Boards where {hero_hand_pretty} wins: {hero_wins}")
    print(f"Boards resulting in a tie: {ties}")
    print(f"Final Equity Equation: ({hero_wins} + {ties} / 2) / {total_boards}")
    print(f"Resulting Equity: {equity:.4%}\n")
    
    # We return the opponent's hand name and our equity against it
    return (villain_cards_str[0][0] + villain_cards_str[1][0], equity)

# --- Main Program ---
# Define the hands for the problem
our_hand = ['As', 'Ac']
opponent_hands = [
    ['Qh', 'Jh'],
    ['Qh', 'Th'],
    ['Qh', '9h']
]

print("Calculating equity for black Aces vs. three red suited hands...")
print("The hand we least want to face is the one that minimizes our equity.\n")

# Calculate equity for each matchup
results = []
for opponent_hand in opponent_hands:
    results.append(calculate_equity_and_print(our_hand, opponent_hand))

# Find the hand that gives us the lowest equity
worst_hand_for_us, min_equity = min(results, key=lambda item: item[1])

print("--- Conclusion ---")
print(f"Comparing the results, the lowest equity for our Aces ({min_equity:.4%}) is against {worst_hand_for_us} suited.")
print(f"Therefore, {worst_hand_for_us} suited is the hand we least like to see.")

# Map the result to the answer choices
if worst_hand_for_us == "QJ":
    final_answer = "A"
elif worst_hand_for_us == "QT":
    final_answer = "B"
elif worst_hand_for_us == "Q9":
    final_answer = "C"
else:
    # This should not be reached based on poker theory and the calculation
    final_answer = "E"

print(f"<<<{final_answer}>>>")