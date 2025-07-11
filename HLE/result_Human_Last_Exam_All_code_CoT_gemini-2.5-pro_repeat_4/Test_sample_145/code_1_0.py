import random
from collections import Counter

def get_hand_rank(hand):
    """
    Evaluates a 7-card poker hand and returns its rank as a comparable tuple.
    A higher tuple value represents a better hand.
    E.g., (8, 14) for Royal Flush, (2, 14, 13, 11) for Two Pair (Aces and Kings, Jack kicker).
    Hand Ranks: 8=SF, 7=Quads, 6=FH, 5=Flush, 4=Straight, 3=Trips, 2=2P, 1=1P, 0=HC
    """
    ranks = sorted([card[0] for card in hand], reverse=True)
    suits = [card[1] for card in hand]
    
    # Check for Flush
    suit_counts = Counter(suits)
    flush_suit = None
    for s, count in suit_counts.items():
        if count >= 5:
            flush_suit = s
            break

    # Check for Straight Flush
    if flush_suit:
        flush_ranks = sorted([r for r, s in hand if s == flush_suit], reverse=True)
        # Check for wheel straight flush (A-5)
        if 14 in flush_ranks and 2 in flush_ranks and 3 in flush_ranks and 4 in flush_ranks and 5 in flush_ranks:
            return (8, 5)  # Rank 8, high card 5
        # Check for regular straight flush
        unique_flush_ranks = sorted(list(set(flush_ranks)), reverse=True)
        for i in range(len(unique_flush_ranks) - 4):
            if unique_flush_ranks[i] - unique_flush_ranks[i+4] == 4:
                return (8, unique_flush_ranks[i])

    # Check for Four of a Kind, Full House, etc.
    rank_counts = Counter(ranks)
    sorted_by_count = sorted(rank_counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
    
    # Four of a Kind
    if sorted_by_count[0][1] == 4:
        quad_rank = sorted_by_count[0][0]
        kicker = max([r for r in ranks if r != quad_rank])
        return (7, quad_rank, kicker)

    # Full House
    if sorted_by_count[0][1] == 3 and sorted_by_count[1][1] >= 2:
        trip_rank = sorted_by_count[0][0]
        pair_rank = max([item[0] for item in sorted_by_count[1:] if item[1] >= 2])
        return (6, trip_rank, pair_rank)
        
    # Return Flush if found
    if flush_suit:
        flush_ranks = sorted([r for r, s in hand if s == flush_suit], reverse=True)
        return (5, tuple(flush_ranks[:5]))

    # Check for Straight
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    # Check for wheel straight (A-5)
    if 14 in unique_ranks and 2 in unique_ranks and 3 in unique_ranks and 4 in unique_ranks and 5 in unique_ranks:
        return (4, 5)
    # Check for regular straight
    for i in range(len(unique_ranks) - 4):
        if unique_ranks[i] - unique_ranks[i+4] == 4:
            return (4, unique_ranks[i])
    
    # Three of a Kind
    if sorted_by_count[0][1] == 3:
        trip_rank = sorted_by_count[0][0]
        kickers = [r for r in ranks if r != trip_rank]
        return (3, trip_rank, kickers[0], kickers[1])
        
    # Two Pair
    if sorted_by_count[0][1] == 2 and sorted_by_count[1][1] == 2:
        pair1, pair2 = sorted_by_count[0][0], sorted_by_count[1][0]
        kicker = max([r for r in ranks if r != pair1 and r != pair2])
        return (2, pair1, pair2, kicker)

    # One Pair
    if sorted_by_count[0][1] == 2:
        pair_rank = sorted_by_count[0][0]
        kickers = [r for r in ranks if r != pair_rank]
        return (1, pair_rank, kickers[0], kickers[1], kickers[2])
        
    # High Card
    return (0, tuple(ranks[:5]))

def calculate_equity(hero_hand, villain_hand, num_simulations):
    """Runs a Monte Carlo simulation to calculate hero's equity."""
    # Create a standard 52-card deck: (rank, suit)
    # Ranks: 2-14 (A=14), Suits: s,c,h,d
    ranks = range(2, 15)
    suits = ['s', 'c', 'h', 'd']
    deck = [(r, s) for r in ranks for s in suits]

    # Remove known cards from the deck
    for card in hero_hand + villain_hand:
        deck.remove(card)

    hero_wins = 0
    villain_wins = 0
    ties = 0

    for _ in range(num_simulations):
        # Draw 5 random cards for the board
        board = random.sample(deck, 5)
        
        # Evaluate hands
        hero_full_hand = hero_hand + board
        villain_full_hand = villain_hand + board
        
        hero_rank = get_hand_rank(hero_full_hand)
        villain_rank = get_hand_rank(villain_full_hand)

        # Compare and tally results
        if hero_rank > villain_rank:
            hero_wins += 1
        elif villain_rank > hero_rank:
            villain_wins += 1
        else:
            ties += 1
            
    hero_equity = (hero_wins + ties / 2) / num_simulations
    return hero_wins, ties, num_simulations, hero_equity

def main():
    """Main function to run the simulations and print results."""
    # Define hands: (Rank, Suit) - A=14, K=13, Q=12, J=11, T=10
    hero_hand = [(14, 's'), (14, 'c')] # Two black Aces
    villain_hands = {
        "QJ suited": [(12, 'h'), (11, 'h')], # Red QJ suited
        "QT suited": [(12, 'h'), (10, 'h')], # Red QT suited
        "Q9 suited": [(12, 'h'), (9, 'h')]   # Red Q9 suited
    }
    num_simulations = 100000
    
    results = {}
    print(f"Running {num_simulations} simulations for each matchup...\n")
    
    for name, v_hand in villain_hands.items():
        wins, ties, total, equity = calculate_equity(hero_hand, v_hand, num_simulations)
        results[name] = equity
        print(f"--- Matchup: Black Aces vs. Red {name} ---")
        print(f"Hero Wins: {wins}, Ties: {ties}, Total Simulations: {total}")
        print(f"Final Equity Equation: ({wins} + {ties}/2) / {total}")
        print(f"Aces' Equity: {equity:.2%}\n")

    # Determine which hand the hero least likes to see
    min_equity_hand = min(results, key=results.get)
    
    print("--- Conclusion ---")
    print("You least like to see the hand that minimizes your equity.")
    print(f"The hand that gives your Aces the lowest equity is: {min_equity_hand}")

if __name__ == '__main__':
    main()
<<<A>>>