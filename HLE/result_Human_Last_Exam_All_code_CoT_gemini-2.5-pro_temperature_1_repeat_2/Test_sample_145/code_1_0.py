import itertools
import random
from collections import Counter

def get_hand_rank(hand):
    """
    Evaluates a 5-card hand and returns a numerical rank.
    Hand is a list of 5 (rank, suit) tuples.
    Rank: 2-14 (Ace=14).
    Returns a tuple, e.g., (8, 14) for a Royal Flush. Higher is better.
    """
    ranks = sorted([card[0] for card in hand], reverse=True)
    suits = {card[1] for card in hand}

    is_flush = len(suits) == 1
    # Check for A-5 low straight (ranks are [14, 5, 4, 3, 2])
    is_straight = (len(set(ranks)) == 5 and ranks[0] - ranks[4] == 4) or (ranks == [14, 5, 4, 3, 2])

    # Straight flush
    if is_straight and is_flush:
        # For A-5 straight flush, the high card is 5
        return (8, 5) if ranks[0] == 14 and ranks[1] == 5 else (8, ranks[0])

    counts = Counter(ranks)
    sorted_counts = sorted(counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
    main_ranks = [item[0] for item in sorted_counts]

    # Four of a kind
    if sorted_counts[0][1] == 4:
        return (7, main_ranks[0], main_ranks[1])
    # Full house
    if sorted_counts[0][1] == 3 and sorted_counts[1][1] == 2:
        return (6, main_ranks[0], main_ranks[1])
    # Flush
    if is_flush:
        return (5, tuple(ranks))
    # Straight
    if is_straight:
        return (4, 5) if ranks[0] == 14 and ranks[1] == 5 else (4, ranks[0])
    # Three of a kind
    if sorted_counts[0][1] == 3:
        return (3, main_ranks[0], main_ranks[1], main_ranks[2])
    # Two pair
    if sorted_counts[0][1] == 2 and sorted_counts[1][1] == 2:
        return (2, main_ranks[0], main_ranks[1], main_ranks[2])
    # One pair
    if sorted_counts[0][1] == 2:
        return (1, main_ranks[0], main_ranks[1], main_ranks[2], main_ranks[3])
    # High card
    return (0, tuple(ranks))

def get_best_hand(seven_cards):
    """Finds the best 5-card hand from a set of 7 cards."""
    return max(get_hand_rank(combo) for combo in itertools.combinations(seven_cards, 5))

def run_simulation(hero_hand, villain_hand, num_simulations=100000):
    """Runs a Monte Carlo simulation to calculate equity."""
    deck = [(r, s) for r in range(2, 15) for s in ['s', 'c', 'h', 'd']]
    
    # Remove known cards
    for card in hero_hand + villain_hand:
        deck.remove(card)

    wins = 0
    ties = 0

    for _ in range(num_simulations):
        random.shuffle(deck)
        board = deck[:5]
        
        hero_best = get_best_hand(hero_hand + board)
        villain_best = get_best_hand(villain_hand + board)

        if hero_best > villain_best:
            wins += 1
        elif hero_best == villain_best:
            ties += 1
            
    return wins, ties, num_simulations

def main():
    """Main function to run the analysis."""
    # Define hands
    # Ranks: J=11, Q=12, K=13, A=14
    hero_hand = [(14, 's'), (14, 'c')] # Black Aces

    villain_hands = {
        "QJ": [(12, 'h'), (11, 'h')], # Red Suited QJ
        "QT": [(12, 'h'), (10, 'h')], # Red Suited QT
        "Q9": [(12, 'h'), (9, 'h')],  # Red Suited Q9
    }

    results = {}
    print("Running simulations to calculate Aces' equity against different hands...\n")

    for name, hand in villain_hands.items():
        wins, ties, total = run_simulation(hero_hand, hand)
        
        win_rate = wins / total
        tie_rate = ties / total
        equity = win_rate + (tie_rate / 2)
        results[name] = equity
        
        print(f"--- Matchup: Black Aces vs Red Suited {name} ---")
        print(f"Final Equation based on {total} simulations:")
        # Print each number in the final equation
        print(f"Equity = (Wins / Total) + (Ties / Total) / 2")
        print(f"Equity = ({wins} / {total}) + ({ties} / {total}) / 2")
        print(f"Equity = {win_rate:.4f} + {tie_rate:.4f} / 2")
        print(f"Resulting Aces' Equity: {equity:.2%}\n")
    
    # Determine the hand that minimizes the Aces' equity
    min_equity_hand = min(results, key=results.get)
    min_equity_value = results[min_equity_hand]

    print("--- Conclusion ---")
    print(f"The hand that minimizes the equity for the black aces is Red Suited {min_equity_hand}, with an equity of {min_equity_value:.2%}.")
    print("Therefore, QJ suited is the hand you least like to see.")

if __name__ == "__main__":
    main()
<<<A>>>