import sys

def solve_poker_equity():
    """
    This script calculates the preflop equity for a pair of black aces
    against three different red suited hands (QJs, QTs, Q9s) in Texas Hold'em.
    It determines which of the opponent's hands is the most difficult for the aces,
    meaning it gives the aces the lowest win probability.

    The script uses the 'treys' library to enumerate all possible boards
    and evaluate the hands efficiently. If you don't have it installed, please run:
    pip install treys
    """
    try:
        from treys import Card, Deck, Evaluator
    except ImportError:
        print("Error: The 'treys' library is required to run this script.", file=sys.stderr)
        print("Please install it by running: pip install treys", file=sys.stderr)
        sys.exit(1)

    # Our hand: Two black aces (Ace of Spades, Ace of Clubs)
    hero_hand_str = ['As', 'Ac']

    # The three opponent hands to test (we'll use diamonds as the red suit)
    villain_hands_options = {
        "QJ": ['Qd', 'Jd'],
        "QT": ['Qd', 'Td'],
        "Q9": ['Qd', '9d']
    }

    print(f"Analyzing matchups for Hero: {hero_hand_str[0]}{hero_hand_str[1]}")
    print("This may take a minute or two to run as it checks every possible outcome...")
    print("-" * 45)

    results = {}
    is_first_calculation = True

    for name, villain_hand_str in villain_hands_options.items():
        # --- Simulation Setup ---
        hero_cards = [Card.new(c) for c in hero_hand_str]
        villain_cards = [Card.new(c) for c in villain_hand_str]

        deck = Deck()
        # Remove the four cards in play from the deck
        for card in hero_cards + villain_cards:
            deck.cards.remove(card)
        
        evaluator = Evaluator()

        # --- Counters ---
        hero_wins = 0
        ties = 0
        total_boards = 0

        # --- Enumerate all possible boards ---
        # The enumerator iterates through all C(48, 5) = 1,712,304 possible boards
        for board in deck.enumerate_boards():
            total_boards += 1
            
            # Evaluate hand strengths. Lower score is better in 'treys'.
            hero_score = evaluator.evaluate(hero_cards, board)
            villain_score = evaluator.evaluate(villain_cards, board)

            if hero_score < villain_score:
                hero_wins += 1
            elif hero_score == villain_score:
                ties += 1
        
        # --- Calculate Equity ---
        equity = (hero_wins + (ties / 2)) / total_boards
        results[name] = equity
        
        print(f"Matchup: {hero_hand_str[0]}{hero_hand_str[1]} vs {villain_hand_str[0]}{villain_hand_str[1]}")

        # As requested, showing the final equation with numbers
        if is_first_calculation:
            print("\nEquity is calculated as: (Hero Wins + (Ties / 2)) / Total Possible Boards")
            print("Final equation for this matchup:")
            print(f"Equity = ({hero_wins} + ({ties} / 2)) / {total_boards}")
            is_first_calculation = False
        
        print(f"Hero's Equity: {equity:.4%}\n")


    # --- Conclusion ---
    min_equity_hand_name = min(results, key=results.get)
    
    print("-" * 45)
    print("Conclusion:")
    print("The equity for our Aces against each hand is:")
    # Sort results by equity value for clearer presentation
    for name, equity in sorted(results.items(), key=lambda item: item[1]):
        print(f"  - vs {name} suited: {equity:.4%}")

    print(f"\nWe least like to see {min_equity_hand_name} suited, because it minimizes our equity.")


if __name__ == '__main__':
    solve_poker_equity()