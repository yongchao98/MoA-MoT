import collections
import itertools

def calculate_equity():
    """
    Calculates the equity of a hero's hand against several villain hands by enumerating all possible boards.
    """

    # --- Hand Evaluation Logic ---
    # Hand rank constants for comparison
    HIGH_CARD = 0
    PAIR = 1
    TWO_PAIR = 2
    THREE_OF_A_KIND = 3
    STRAIGHT = 4
    FLUSH = 5
    FULL_HOUSE = 6
    FOUR_OF_A_KIND = 7
    STRAIGHT_FLUSH = 8

    # Card rank mapping (T=10, J=11, Q=12, K=13, A=14)
    RANK_MAP = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}

    def evaluate_hand(cards):
        """
        Evaluates a 7-card hand and returns its strength as a comparable tuple.
        Example card format: 'As', 'Td', '7c'.
        Return format: (hand_type, kicker1, kicker2, ...)
        """
        ranks = sorted([RANK_MAP[c[0]] for c in cards], reverse=True)
        suits = [c[1] for c in cards]

        # Check for flush and straight flush
        suit_counts = collections.Counter(suits)
        flush_suit = next((s for s, count in suit_counts.items() if count >= 5), None)

        is_flush = flush_suit is not None
        is_straight = False
        straight_high_card = -1

        eval_ranks = ranks
        if is_flush:
            flush_ranks = sorted([RANK_MAP[c[0]] for c in cards if c[1] == flush_suit], reverse=True)
            eval_ranks = flush_ranks # Check for straight within the flush

        unique_ranks = sorted(list(set(eval_ranks)), reverse=True)
        if len(unique_ranks) >= 5:
            for i in range(len(unique_ranks) - 4):
                if unique_ranks[i] - unique_ranks[i+4] == 4:
                    is_straight = True
                    straight_high_card = unique_ranks[i]
                    break
            # Ace-low straight check (A,2,3,4,5)
            if not is_straight and set([14, 2, 3, 4, 5]).issubset(set(unique_ranks)):
                is_straight = True
                straight_high_card = 5

        if is_straight and is_flush:
            return (STRAIGHT_FLUSH, straight_high_card)
        if is_flush:
            return (FLUSH, tuple(flush_ranks[:5]))
        if is_straight:
            return (STRAIGHT, straight_high_card)

        # Check for pairs, trips, quads, etc.
        rank_counts = collections.Counter(ranks)
        # Sort by count first, then rank to handle full houses correctly
        counts = sorted(rank_counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
        
        # Four of a Kind
        if counts[0][1] == 4:
            quad_rank = counts[0][0]
            kicker = next(r for r in ranks if r != quad_rank)
            return (FOUR_OF_A_KIND, quad_rank, kicker)

        # Full House
        if counts[0][1] == 3 and counts[1][1] >= 2:
            return (FULL_HOUSE, counts[0][0], counts[1][0])

        # Three of a Kind
        if counts[0][1] == 3:
            trip_rank = counts[0][0]
            kickers = [r for r in ranks if r != trip_rank]
            return (THREE_OF_A_KIND, trip_rank, tuple(kickers[:2]))

        # Two Pair
        if counts[0][1] == 2 and len(counts) > 1 and counts[1][1] == 2:
            p1, p2 = counts[0][0], counts[1][0]
            kicker = next(r for r in ranks if r != p1 and r != p2)
            return (TWO_PAIR, p1, p2, kicker)

        # Pair
        if counts[0][1] == 2:
            pair_rank = counts[0][0]
            kickers = [r for r in ranks if r != pair_rank]
            return (PAIR, pair_rank, tuple(kickers[:3]))

        # High Card
        return (HIGH_CARD, tuple(ranks[:5]))

    # --- Simulation Logic ---
    hero_hand = ['As', 'Ac']
    villain_hands = {
        "QJ suited": ['Qh', 'Jh'],
        "QT suited": ['Qh', 'Th'],
        "Q9 suited": ['Qh', '9h']
    }
    
    ranks_str = '23456789TJQKA'
    suits_str = 'shdc'
    full_deck = {r + s for r in ranks_str for s in suits_str}

    results = {}
    print("Calculating equities for Hero (AsAc)... This may take a minute.")

    for name, villain_hand in villain_hands.items():
        deck = list(full_deck - set(hero_hand) - set(villain_hand))
        
        wins, ties, total_boards = 0, 0, 0
        
        for board in itertools.combinations(deck, 5):
            total_boards += 1
            hero_full_hand = hero_hand + list(board)
            villain_full_hand = villain_hand + list(board)
            
            hero_score = evaluate_hand(hero_full_hand)
            villain_score = evaluate_hand(villain_full_hand)

            if hero_score > villain_score:
                wins += 1
            elif hero_score == villain_score:
                ties += 1

        equity = (wins + ties / 2) / total_boards
        results[name] = equity

    # --- Output Results ---
    print("\n--- Final Equity Results ---")
    min_equity_hand = None
    min_equity_value = 1.0

    for hand, equity in results.items():
        print(f"Equity for AsAc vs {hand}: {equity:.2%}")
        if equity < min_equity_value:
            min_equity_value = equity
            min_equity_hand = hand

    print(f"\nThe hand you least like to see is {min_equity_hand}.")
    print("This is because it gives your aces the lowest probability of winning.")
    
    return min_equity_hand


# Run the calculation and determine the final answer choice
final_hand = calculate_equity()
if "QJ" in final_hand:
    answer = "A"
elif "QT" in final_hand:
    answer = "B"
elif "Q9" in final_hand:
    answer = "C"
else:
    answer = "E" # Should not happen

print(f"\nFinal Answer Choice: {answer}")
print(f'<<<{answer}>>>')
