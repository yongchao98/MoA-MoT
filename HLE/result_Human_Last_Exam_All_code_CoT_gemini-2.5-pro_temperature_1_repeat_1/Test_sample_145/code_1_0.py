import random
from collections import Counter
import itertools

def get_hand_rank(hand):
    """
    Evaluates a 7-card hand and returns its rank as a comparable tuple.
    A higher tuple value indicates a better hand.
    Hand Ranks: 9=StraightFlush, 8=Quads, 7=FullHouse, 6=Flush, 5=Straight,
                  4=Trips, 3=TwoPair, 2=Pair, 1=HighCard.
    """
    vals = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}
    
    ranks = sorted([vals[c[0]] for c in hand], reverse=True)
    suits = [c[1] for c in hand]
    
    # Check for Flush
    suit_counts = Counter(suits)
    flush_suit = next((s for s, count in suit_counts.items() if count >= 5), None)
    
    is_flush = flush_suit is not None
    
    # Check for Straight
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    is_straight = False
    straight_high_card = -1
    
    # Ace-low straight check
    if set([14, 2, 3, 4, 5]).issubset(set(unique_ranks)):
        is_straight = True
        straight_high_card = 5
        
    # Regular straight check
    for i in range(len(unique_ranks) - 4):
        if unique_ranks[i] - unique_ranks[i+4] == 4:
            is_straight = True
            straight_high_card = max(straight_high_card, unique_ranks[i])
            break

    # Evaluate Straight Flush
    if is_flush and is_straight:
        flush_ranks = sorted([vals[c[0]] for c in hand if c[1] == flush_suit], reverse=True)
        # Check if the flush cards form a straight
        if set([14, 2, 3, 4, 5]).issubset(set(flush_ranks)):
            return (9, 5) # Straight flush (A-5)
        for i in range(len(flush_ranks) - 4):
            if flush_ranks[i] - flush_ranks[i+4] == 4:
                return (9, flush_ranks[i]) # Straight flush

    # --- Check for Quads, Full House, Trips, Pairs ---
    rank_counts = Counter(ranks)
    # Sort ranks by count, then by rank value
    sorted_ranks_by_count = sorted(rank_counts.keys(), key=lambda k: (rank_counts[k], k), reverse=True)

    counts = sorted(rank_counts.values(), reverse=True)
    
    # Four of a kind
    if counts[0] == 4:
        # sorted_ranks_by_count[0] is the rank of the quads
        # sorted_ranks_by_count[1] is the highest kicker
        return (8, sorted_ranks_by_count[0], sorted_ranks_by_count[1])
    
    # Full House
    if counts[0] == 3 and counts[1] >= 2:
        # sorted_ranks_by_count[0] is the trips, [1] is the pair
        return (7, sorted_ranks_by_count[0], sorted_ranks_by_count[1])

    # Flush (if not a straight flush)
    if is_flush:
        flush_ranks = sorted([vals[c[0]] for c in hand if c[1] == flush_suit], reverse=True)
        return (6, tuple(flush_ranks[:5]))

    # Straight (if not a flush)
    if is_straight:
        return (5, straight_high_card)

    # Three of a kind
    if counts[0] == 3:
        # [0] is trips, [1] and [2] are kickers
        return (4, sorted_ranks_by_count[0], sorted_ranks_by_count[1], sorted_ranks_by_count[2])

    # Two Pair
    if counts[0] == 2 and counts[1] == 2:
        # [0] and [1] are pairs, [2] is kicker
        return (3, sorted_ranks_by_count[0], sorted_ranks_by_count[1], sorted_ranks_by_count[2])

    # One Pair
    if counts[0] == 2:
        # [0] is pair, [1, 2, 3] are kickers
        return (2, sorted_ranks_by_count[0], sorted_ranks_by_count[1], sorted_ranks_by_count[2], sorted_ranks_by_count[3])
    
    # High Card
    return (1, tuple(sorted_ranks_by_count[:5]))


def calculate_equity(hero_hand, villain_hand, simulations=100000):
    """Calculates hero's equity via Monte Carlo simulation."""
    deck = {r + s for r in '23456789TJQKA' for s in 'shdc'}
    deck -= set(hero_hand)
    deck -= set(villain_hand)
    deck_list = list(deck)
    
    wins, ties = 0, 0
    
    for _ in range(simulations):
        board = random.sample(deck_list, 5)
        
        hero_rank = get_hand_rank(hero_hand + tuple(board))
        villain_rank = get_hand_rank(villain_hand + tuple(board))

        if hero_rank > villain_rank:
            wins += 1
        elif hero_rank == villain_rank:
            ties += 1

    equity = (wins + ties / 2) / simulations
    return equity


def main():
    """Main function to run the simulation and print results."""
    hero_hand = ('As', 'Ac') # Two black aces

    # Opponent hands are red suited (we'll use hearts for simulation)
    villain_hands = {
        "QJ suited": ('Qh', 'Jh'),
        "QT suited": ('Qh', 'Th'),
        "Q9 suited": ('Qh', '9h')
    }
    
    print("Calculating equity for Aces against different red suited hands.")
    print("This may take a moment...\n")

    results = {}
    for name, hand in villain_hands.items():
        equity = calculate_equity(hero_hand, hand)
        results[name] = equity
        # This is the "equation" showing each calculated value
        print(f"Equity for {hero_hand[0]}{hero_hand[1]} vs {hand[0]}{hand[1]}: {equity:.2%}")

    # Find the hand that minimizes the hero's equity
    min_equity_hand = min(results, key=results.get)
    
    print(f"\nThe hand you least like to see is the one that minimizes your equity.")
    print(f"The hand that gives the Aces the lowest equity is: {min_equity_hand}")


if __name__ == "__main__":
    main()

<<<A>>>