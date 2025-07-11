import itertools
from deuces import Card, Deck, Evaluator

def run_equity_simulation(hero_hand_str, villain_hand_str, num_simulations):
    """
    Runs a Monte Carlo simulation to calculate poker equity.
    """
    hero_hand = [Card.new(c) for c in hero_hand_str]
    villain_hand = [Card.new(c) for c in villain_hand_str]
    
    # These are the cards we know are not in the deck for the board
    dead_cards = hero_hand + villain_hand

    wins = 0
    ties = 0
    losses = 0
    
    # Create the evaluator instance
    evaluator = Evaluator()

    # Create a deck and remove the known cards once
    deck = Deck()
    for card in dead_cards:
        deck.cards.remove(card)
    
    # Get the list of remaining cards for combinations
    initial_deck_list = list(deck.cards)

    # Use itertools.combinations for a full enumeration if num_simulations is large enough
    # or a Monte Carlo simulation otherwise. For speed, we will do Monte Carlo.
    # Note: Full enumeration would be `itertools.combinations(initial_deck_list, 5)`
    
    for i in range(num_simulations):
        # Draw a 5-card board randomly
        deck.shuffle()
        board = deck.draw(5)

        # Evaluate the hero's hand and villain's hand
        # Lower score is better in the deuces library
        hero_score = evaluator.evaluate(board, hero_hand)
        villain_score = evaluator.evaluate(board, villain_hand)

        if hero_score < villain_score:
            wins += 1
        elif hero_score > villain_score:
            losses += 1
        else:
            ties += 1
            
    # Calculate equity
    equity = (wins + 0.5 * ties) / num_simulations
    
    # Print the detailed results including the equation
    print(f"Matchup: {' '.join(hero_hand_str)} vs {' '.join(villain_hand_str)}")
    print(f"Simulated {num_simulations} hands.")
    print(f"Hero Wins: {wins}, Ties: {ties}, Losses: {losses}")
    print(f"Hero Equity = ({wins} + 0.5 * {ties}) / {num_simulations} = {equity:.4f}")
    
    return equity

if __name__ == '__main__':
    # Hero's hand: Two black aces
    hero_hand = ["As", "Ac"]
    
    # Villain hands to test (red suited)
    villain_hands = {
        "QJ": ["Qh", "Jh"],
        "QT": ["Qh", "Th"],
        "Q9": ["Qh", "9h"],
    }
    
    # Number of simulations for each matchup
    num_simulations = 200000
    
    results = {}
    
    print("Starting poker equity simulations...")
    print("-" * 40)
    
    for name, hand_str in villain_hands.items():
        equity = run_equity_simulation(hero_hand, hand_str, num_simulations)
        results[name] = equity
        print("-" * 40)
        
    # Determine which hand is the worst for the hero (lowest equity)
    worst_hand_for_hero = min(results, key=results.get)
    min_equity = results[worst_hand_for_hero]
    
    print("\n--- Conclusion ---")
    print(f"The equity for {hero_hand[0]}{hero_hand[1]} is lowest against {villain_hands[worst_hand_for_hero][0]}{villain_hands[worst_hand_for_hero][1]} ({worst_hand_for_hero} suited).")
    print(f"This is because connected cards (like QJ) have a higher potential to make straights and straight flushes than gapped cards (like QT or Q9).")
    print(f"Therefore, the hand you least like to see is {worst_hand_for_hero} suited.")

    # Map the result to the answer choice
    answer_map = {"QJ": "A", "QT": "B", "Q9": "C"}
    final_answer = answer_map.get(worst_hand_for_hero, "Error")
    
    print(f"\nFinal Answer Choice: {final_answer}")
    print(f"<<<{final_answer}>>>")