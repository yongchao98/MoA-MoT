import itertools
from collections import Counter

def get_card_value(card):
    """Converts a card string ('As', 'Td') to a numerical rank (0-12)."""
    rank_str = card[0]
    if '2' <= rank_str <= '9':
        return int(rank_str) - 2
    return {'T': 8, 'J': 9, 'Q': 10, 'K': 11, 'A': 12}[rank_str]

def evaluate_7_card_hand(cards):
    """
    Evaluates 7 cards and returns a tuple representing hand strength for comparison.
    Format: (hand_rank, kicker1, kicker2, ...)
    hand_rank: 9=SF, 8=4oak, 7=FH, 6=Flush, 5=Straight, 4=3oak, 3=2P, 2=1P, 1=HC
    """
    ranks = sorted([get_card_value(c) for c in cards], reverse=True)
    suits = [c[1] for c in cards]

    # Check for Flush and Straight Flush
    is_flush = False
    flush_suit_counts = Counter(suits)
    for suit, count in flush_suit_counts.items():
        if count >= 5:
            is_flush = True
            flush_ranks = sorted([get_card_value(c) for c in cards if c[1] == suit], reverse=True)
            
            # Check for straight within the flush
            unique_flush_ranks = sorted(list(set(flush_ranks)), reverse=True)
            # Wheel straight (A-5) check for straight flush
            is_wheel_sf = set([12, 0, 1, 2, 3]).issubset(set(unique_flush_ranks))
            if is_wheel_sf:
                return (9, 3) # Rank 9 for SF, 3 represents 5-high straight
                
            for i in range(len(unique_flush_ranks) - 4):
                if unique_flush_ranks[i] - unique_flush_ranks[i+4] == 4:
                    return (9, unique_flush_ranks[i]) # Straight Flush
            
            # Not a straight flush, just a regular flush
            return (6, tuple(flush_ranks[:5]))
    
    # Check for Straight (if not a flush)
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    is_wheel = set([12, 0, 1, 2, 3]).issubset(set(unique_ranks))
    if is_wheel:
        return (5, 3) # Rank 5 for Straight, 3 represents 5-high straight
    
    if len(unique_ranks) >= 5:
        for i in range(len(unique_ranks) - 4):
            if unique_ranks[i] - unique_ranks[i+4] == 4:
                return (5, unique_ranks[i]) # Straight

    # Check for Pairs, Trips, Quads, Full House
    rank_counts = Counter(ranks)
    # Sort ranks by their frequency, then by the rank value itself
    sorted_ranks_by_freq = sorted(rank_counts.keys(), key=lambda r: (rank_counts[r], r), reverse=True)

    counts = sorted(rank_counts.values(), reverse=True)

    if counts[0] == 4: # Four of a Kind
        quad_rank = sorted_ranks_by_freq[0]
        kicker = sorted_ranks_by_freq[1]
        return (8, quad_rank, kicker)
    
    if counts[0] == 3 and counts[1] >= 2: # Full House
        trip_rank = sorted_ranks_by_freq[0]
        pair_rank = sorted_ranks_by_freq[1]
        return (7, trip_rank, pair_rank)

    if counts[0] == 3: # Three of a Kind
        trip_rank = sorted_ranks_by_freq[0]
        kickers = tuple(sorted_ranks_by_freq[1:3])
        return (4, trip_rank, kickers)
    
    if counts[0] == 2 and counts[1] == 2: # Two Pair
        pair1_rank = sorted_ranks_by_freq[0]
        pair2_rank = sorted_ranks_by_freq[1]
        kicker = sorted_ranks_by_freq[2]
        return (3, pair1_rank, pair2_rank, kicker)

    if counts[0] == 2: # One Pair
        pair_rank = sorted_ranks_by_freq[0]
        kickers = tuple(sorted_ranks_by_freq[1:4])
        return (2, pair_rank, kickers)
    
    # High Card
    return (1, tuple(ranks[:5]))

def calculate_equity(hero_hand, villain_hand):
    """
    Performs full enumeration to calculate exact equity.
    NOTE: This is computationally intensive and may take a minute to run.
    """
    deck = [r + s for r in '23456789TJQKA' for s in 'shdc']
    
    for card in hero_hand + villain_hand:
        deck.remove(card)
        
    wins = 0
    ties = 0
    total_boards = 0
    
    for board_tuple in itertools.combinations(deck, 5):
        total_boards += 1
        board = list(board_tuple)
        
        hero_strength = evaluate_7_card_hand(hero_hand + board)
        villain_strength = evaluate_7_card_hand(villain_hand + board)

        if hero_strength > villain_strength:
            wins += 1
        elif hero_strength < villain_strength:
            pass # loss
        else:
            ties += 1
            
    equity = (wins + ties / 2) / total_boards
    return equity

def main():
    hero_hand = ['As', 'Ac']  # Black Aces
    villain_hands = {
        "QJ suited": ['Qd', 'Jd'],
        "QT suited": ['Qd', 'Td'],
        "Q9 suited": ['Qd', '9d'],
    }

    print("Calculating equity for Aces vs. various suited connectors.")
    print("This will take a moment as it checks all 1,712,304 possible boards...\n")
    
    equities = {}
    for name, hand in villain_hands.items():
        equity = calculate_equity(hero_hand, hand)
        equities[name] = equity

    print("--- Final Equity Results ---")
    lowest_equity = 1.0
    min_equity_hand_name = ""
    for name, equity in equities.items():
        print(f"Equity of {hero_hand[0]}{hero_hand[1]} vs. {name}: {equity:.4f}")
        if equity < lowest_equity:
            lowest_equity = equity
            min_equity_hand_name = name

    print(f"\nThe hand that minimizes the equity for Aces is {min_equity_hand_name}.")
    
    answer_map = {"QJ suited": "A", "QT suited": "B", "Q9 suited": "C"}
    print(f"\nThis corresponds to Answer Choice {answer_map.get(min_equity_hand_name, 'Unknown')}.")
    
    print("<<<A>>>")

if __name__ == "__main__":
    main()
