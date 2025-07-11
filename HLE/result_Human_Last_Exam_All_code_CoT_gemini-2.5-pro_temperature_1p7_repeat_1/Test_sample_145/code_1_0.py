import collections
from itertools import combinations

def calculate_equity():
    """
    Calculates the equity of Aces against three different suited hands (QJs, QTs, Q9s)
    and determines which opponent hand is the most dangerous.
    """

    # --- Hand Evaluation Logic ---
    RANKS = '23456789TJQKA'
    RANK_MAP = {rank: i for i, rank in enumerate(RANKS)}

    def evaluate_hand(cards):
        """
        Evaluates a 7-card hand and returns a tuple representing its strength.
        The tuple can be directly compared to other hand strength tuples.
        Example: (8, 12) for a Royal Flush vs (7, 10, 11) for Four of a Kind (Jacks) with a King kicker.
        """
        ranks__int = sorted([RANK_MAP[c[0]] for c in cards], reverse=True)
        suits = [c[1] for c in cards]
        
        rank_counts = collections.Counter(ranks_int)
        suit_counts = collections.Counter(suits)

        # Check for flush
        is_flush = False
        flush_suit = None
        flush_ranks = []
        for suit, count in suit_counts.items():
            if count >= 5:
                is_flush = True
                flush_suit = suit
                flush_ranks = sorted([RANK_MAP[c[0]] for c in cards if c[1] == flush_suit], reverse=True)
                break
        
        # Check for straight
        unique_ranks = sorted(list(set(ranks_int)), reverse=True)
        is_straight = False
        straight_high_rank = -1
        
        # Ace-low (wheel) straight check: A, 5, 4, 3, 2
        if all(x in unique_ranks for x in [RANK_MAP['A'], RANK_MAP['5'], RANK_MAP['4'], RANK_MAP['3'], RANK_MAP['2']]):
            is_straight = True
            straight_high_rank = RANK_MAP['5']

        for i in range(len(unique_ranks) - 4):
            if not is_straight and unique_ranks[i] - unique_ranks[i+4] == 4:
                is_straight = True
                straight_high_rank = unique_ranks[i]
                break
        
        # Determine hand type and return comparable tuple
        # 1. Straight Flush
        if is_flush and is_straight:
            # Check if the flush cards form a straight
            unique_flush_ranks = sorted(list(set(flush_ranks)), reverse=True)
            # Ace-low straight flush
            if all(x in unique_flush_ranks for x in [RANK_MAP['A'], RANK_MAP['5'], RANK_MAP['4'], RANK_MAP['3'], RANK_MAP['2']]):
                return (8, RANK_MAP['5'])
            # Regular straight flush
            for i in range(len(unique_flush_ranks) - 4):
                if unique_flush_ranks[i] - unique_flush_ranks[i+4] == 4:
                    return (8, unique_flush_ranks[i])

        # 2. Four of a Kind
        if 4 in rank_counts.values():
            quad_rank = [r for r, c in rank_counts.items() if c == 4][0]
            kicker = max([r for r in ranks_int if r != quad_rank])
            return (7, quad_rank, kicker)

        # 3. Full House
        threes = [r for r, c in rank_counts.items() if c == 3]
        twos = [r for r, c in rank_counts.items() if c == 2]
        if len(threes) > 1 or (len(threes) == 1 and len(twos) >= 1):
            trip_rank = max(threes)
            pair_pool = [r for r in threes if r != trip_rank] + twos
            pair_rank = max(pair_pool)
            return (6, trip_rank, pair_rank)

        # 4. Flush
        if is_flush:
            return (5, tuple(flush_ranks[:5]))

        # 5. Straight
        if is_straight:
            return (4, straight_high_rank)

        # 6. Three of a Kind
        if len(threes) == 1:
            trip_rank = threes[0]
            kickers = [r for r in ranks_int if r != trip_rank]
            return (3, trip_rank, kickers[0], kickers[1])
        
        # 7. Two Pair
        if len(twos) >= 2:
            high_pair = max(twos)
            low_pair = max([r for r in twos if r != high_pair])
            kicker = max([r for r in ranks_int if r != high_pair and r != low_pair])
            return (2, high_pair, low_pair, kicker)

        # 8. One Pair
        if len(twos) == 1:
            pair_rank = twos[0]
            kickers = [r for r in ranks_int if r != pair_rank]
            return (1, pair_rank, kickers[0], kickers[1], kickers[2])
            
        # 9. High Card
        return (0, tuple(ranks_int[:5]))


    # --- Main Simulation Logic ---
    hero_hand = ['As', 'Ac']
    villain_hands_map = {
        "QJ suited": ['Qd', 'Jd'],
        "QT suited": ['Qd', 'Td'],
        "Q9 suited": ['Qd', '9d'],
    }
    
    all_ranks = '23456789TJQKA'
    all_suits = 'shdc'
    full_deck = [r + s for r in all_ranks for s in all_suits]
    
    results = {}

    for name, villain_hand in villain_hands_map.items():
        print(f"--- Calculating Equity for {''.join(hero_hand)} vs {''.join(villain_hand)} ---")
        
        known_cards = set(hero_hand + villain_hand)
        remaining_deck = [card for card in full_deck if card not in known_cards]
        
        wins = 0
        ties = 0
        losses = 0
        
        total_boards = 0
        for board in combinations(remaining_deck, 5):
            total_boards += 1
            hero_7_cards = hero_hand + list(board)
            villain_7_cards = villain_hand + list(board)
            
            hero_rank = evaluate_hand(hero_7_cards)
            villain_rank = evaluate_hand(villain_7_cards)
            
            if hero_rank > villain_rank:
                wins += 1
            elif villain_rank > hero_rank:
                losses += 1
            else:
                ties += 1

        equity = (wins + (ties / 2)) / total_boards
        results[name] = equity

        print(f"Total boards analyzed: {total_boards}")
        print(f"Aces Win Count: {wins}")
        print(f"{name} Win Count: {losses}")
        print(f"Tie Count: {ties}")
        print(f"Aces Equity Equation: ({wins} + {ties} / 2) / {total_boards}")
        print(f"Resulting Equity for Aces: {equity:.2%}\n")
    
    min_equity_hand = min(results, key=results.get)
    print("--- Conclusion ---")
    print(f"The hand that minimizes your equity (the one you least like to see) is {min_equity_hand}.")


calculate_equity()