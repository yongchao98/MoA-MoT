import itertools
from collections import Counter

def solve():
    """
    Calculates the equity for Aces against three different suited hands
    in Texas Hold'em to determine which is the most dangerous.
    """
    
    # Map card ranks to numerical values for easy comparison.
    # T=10, J=11, Q=12, K=13, A=14
    val_map = {r: i for i, r in enumerate('23456789TJQKA', 2)}

    def get_hand_rank(hand):
        """
        Evaluates a 7-card hand and returns a tuple representing its strength.
        A higher tuple value means a better hand.
        E.g., Royal Flush > Straight Flush > Four of a Kind, etc.
        """
        ranks = sorted([val_map[r[0]] for r in hand], reverse=True)
        suits = [s[1] for s in hand]

        suit_counts = Counter(suits)
        flush_suit = next((s for s, c in suit_counts.items() if c >= 5), None)

        unique_ranks = sorted(list(set(ranks)), reverse=True)
        
        is_straight = False
        # Check for Ace-low straight (wheel) first
        if set([14, 2, 3, 4, 5]).issubset(unique_ranks):
            straight_high_rank = 5
            is_straight = True
        else:
            # Check for other straights
            for i in range(len(unique_ranks) - 4):
                if unique_ranks[i] - unique_ranks[i+4] == 4:
                    straight_high_rank = unique_ranks[i]
                    is_straight = True
                    break

        # 1. Straight Flush (and Royal Flush)
        if is_straight and flush_suit:
            flush_ranks = sorted([val_map[r[0]] for r in hand if r[1] == flush_suit], reverse=True)
            # Check for Ace-low straight flush
            if set([14, 2, 3, 4, 5]).issubset(flush_ranks):
                return (8, 5) # Rank 8, high card 5
            for i in range(len(flush_ranks) - 4):
                if flush_ranks[i] - flush_ranks[i+4] == 4:
                    return (8, flush_ranks[i])

        rank_counts = Counter(ranks)
        sorted_counts = sorted(rank_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        
        # 2. Four of a Kind
        if sorted_counts[0][1] == 4:
            quad_rank = sorted_counts[0][0]
            kicker = max(r for r in ranks if r != quad_rank)
            return (7, quad_rank, kicker)

        # 3. Full House
        if sorted_counts[0][1] == 3 and sorted_counts[1][1] >= 2:
            return (6, sorted_counts[0][0], sorted_counts[1][0])

        # 4. Flush
        if flush_suit:
            flush_ranks = sorted([val_map[r[0]] for r in hand if r[1] == flush_suit], reverse=True)
            return (5, tuple(flush_ranks[:5]))

        # 5. Straight
        if is_straight:
            return (4, straight_high_rank)

        # 6. Three of a Kind
        if sorted_counts[0][1] == 3:
            trips_rank = sorted_counts[0][0]
            kickers = [r for r in ranks if r != trips_rank]
            return (3, trips_rank, kickers[0], kickers[1])

        # 7. Two Pair
        if sorted_counts[0][1] == 2 and sorted_counts[1][1] == 2:
            p1, p2 = sorted_counts[0][0], sorted_counts[1][0]
            kicker = max(r for r in ranks if r not in [p1, p2])
            return (2, p1, p2, kicker)

        # 8. One Pair
        if sorted_counts[0][1] == 2:
            pair_rank = sorted_counts[0][0]
            kickers = [r for r in ranks if r != pair_rank]
            return (1, pair_rank, kickers[0], kickers[1], kickers[2])

        # 9. High Card
        return (0, tuple(ranks[:5]))

    # --- Main Logic ---
    
    ranks = "23456789TJQKA"
    suits = "shdc" # Spades, Hearts, Diamonds, Clubs
    deck = [r + s for r in ranks for s in suits]

    hero_hand = ('As', 'Ac') # Two black aces

    # The three hands to test against
    villain_hands_to_test = {
        "QJ": ('Qh', 'Jh'),
        "QT": ('Qh', 'Th'),
        "Q9": ('Qh', '9h'),
    }

    results = {}
    
    print("Calculating equities... (This may take up to a minute)")
    
    for name, villain_hand in villain_hands_to_test.items():
        # Create the remaining deck after removing hole cards
        remaining_deck = [card for card in deck if card not in hero_hand and card not in villain_hand]
        
        hero_wins = 0
        villain_wins = 0
        ties = 0

        # Iterate through all possible 5-card boards
        all_boards = itertools.combinations(remaining_deck, 5)
        total_boards = 1712304 # C(48, 5)

        for board in all_boards:
            hero_full_hand = hero_hand + board
            villain_full_hand = villain_hand + board
            
            hero_rank = get_hand_rank(hero_full_hand)
            villain_rank = get_hand_rank(villain_full_hand)

            if hero_rank > villain_rank:
                hero_wins += 1
            elif villain_rank > hero_rank:
                villain_wins += 1
            else:
                ties += 1
        
        # Calculate equity: (wins + ties/2) / total
        hero_equity = (hero_wins + ties / 2) / total_boards * 100
        results[name] = hero_equity
        
        print(f"\nMatchup: Black Aces vs. Red {name}s")
        print(f"Total Boards Analyzed: {total_boards}")
        print(f"Aces' Wins: {hero_wins}")
        print(f"Aces' Losses: {villain_wins}")
        print(f"Ties: {ties}")
        print(f"Aces' Equity = ({hero_wins} + {ties} / 2) / {total_boards} = {hero_equity:.2f}%")

    # Find the hand that minimizes the Aces' equity
    worst_opponent = min(results, key=results.get)
    min_equity = results[worst_opponent]

    print("\n--- Conclusion ---")
    print(f"The hand that minimizes the equity of the two black aces is {worst_opponent} suited.")
    print(f"Against {worst_opponent}s, the aces have their lowest equity of {min_equity:.2f}%.")
    print("Therefore, this is the hand you would least like to see.")

solve()
<<<A>>>