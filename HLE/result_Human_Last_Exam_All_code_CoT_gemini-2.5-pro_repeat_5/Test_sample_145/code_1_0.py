import itertools
try:
    from deuces import Card, Evaluator
except ImportError:
    print("The 'deuces' library is required. Please install it by running: pip install deuces")
    exit()

def calculate_equity(hero_hand_str, villain_hand_str):
    """
    Calculates the equity of hero_hand vs villain_hand by enumerating all possible boards.
    """
    # Initialize the hand evaluator
    evaluator = Evaluator()

    # Create card objects for the hands
    hero_hand = [Card.new(c) for c in hero_hand_str]
    villain_hand = [Card.new(c) for c in villain_hand_str]

    # Create a full deck of cards and remove the known cards (hero + villain)
    dead_cards = hero_hand + villain_hand
    deck = [c for c in Card.GetFullDeck() if c not in dead_cards]

    # Initialize counters for wins, losses, and ties
    hero_wins = 0
    villain_wins = 0
    ties = 0
    total_boards = 0

    # Iterate through all possible 5-card boards from the remaining deck
    for board in itertools.combinations(deck, 5):
        total_boards += 1

        # Evaluate the 7-card hand for both players
        # The evaluator returns a numerical score; lower is better.
        hero_score = evaluator.evaluate(hero_hand, list(board))
        villain_score = evaluator.evaluate(villain_hand, list(board))

        # Compare scores and update counters
        if hero_score < villain_score:
            hero_wins += 1
        elif villain_score < hero_score:
            villain_wins += 1
        else:
            ties += 1

    # Calculate the hero's equity
    equity = (hero_wins + ties / 2) / total_boards

    return {
        "hero_wins": hero_wins,
        "villain_wins": villain_wins,
        "ties": ties,
        "total_boards": total_boards,
        "equity": equity
    }

def main():
    """
    Main function to run the equity calculations and determine the answer.
    """
    # Define the hands for the scenario
    # Hero has two black aces (Spades and Clubs)
    hero_hand_str = ['As', 'Ac']

    # Villain has one of three red suited hands (we'll use Hearts)
    villain_hands = {
        "QJ suited": ['Qh', 'Jh'],
        "QT suited": ['Qh', 'Th'],
        "Q9 suited": ['Qh', '9h']
    }

    results = {}

    print("Calculating preflop equity for Aces vs various red suited hands.")
    print("This may take a minute as it enumerates over 1.7 million possible boards for each matchup.\n")

    # Calculate equity for each matchup
    for name, hand in villain_hands.items():
        print(f"Calculating for AsAc vs {hand[0]}{hand[1]}...")
        result = calculate_equity(hero_hand_str, hand)
        results[name] = result['equity']

        print(f"Result for matchup: Aces vs {name} ({hand[0]}{hand[1]})")
        print(f"  - Total Boards Checked: {result['total_boards']:,}")
        print(f"  - Aces Win Count:    {result['hero_wins']:,}")
        print(f"  - Opponent Win Count:  {result['villain_wins']:,}")
        print(f"  - Tie Count:         {result['ties']:,}")
        print(f"  - Aces Equity Equation: ({result['hero_wins']} + {result['ties']} / 2) / {result['total_boards']}")
        print(f"  - Aces Equity:    {result['equity']:.2%}\n")

    # Find the hand that minimizes the Aces' equity
    min_equity = 1.0
    worst_opponent = None
    for name, equity in results.items():
        if equity < min_equity:
            min_equity = equity
            worst_opponent = name

    print("--- Conclusion ---")
    print("To find the hand we 'least like to see', we look for the one that gives our Aces the lowest equity.")
    print(f"Comparing the final equities for the Aces:")
    for name, equity in sorted(results.items(), key=lambda item: item[1]):
        print(f"  - vs {name}: {equity:.2%}")

    print(f"\nThe hand that minimizes the Aces' equity is {worst_opponent}.")
    print("Therefore, this is the hand you least like to see.")


if __name__ == "__main__":
    main()
<<<A>>>