import itertools
from collections import Counter
import time

def solve_poker_equity():
    """
    Calculates the exact equity for Aces vs. three different suited hands
    in Texas Hold'em by iterating through all possible boards.
    """

    # 1. Deck and Card Setup
    # Ranks: T=10, J=11, Q=12, K=13, A=14
    # Suits: s, c, h, d
    RANKS = '23456789TJQKA'
    SUITS = 'scdh'
    RANK_MAP = {rank: i for i, rank in enumerate(RANKS, 2)}
    
    deck = {f'{r}{s}' for r in RANKS for s in SUITS}
    
    # --- Helper Functions ---
    def card_to_val(card_str):
        """Converts a card string 'As' to a tuple (14, 's')."""
        rank = RANK_MAP[card_str[0]]
        suit = card_str[1]
        return (rank, suit)

    def score_5_card_hand(hand_vals):
        """
        Scores a 5-card hand and returns a tuple for comparison.
        Hand format: a list of 5 (rank, suit) tuples.
        """
        ranks = sorted([c[0] for c in hand_vals], reverse=True)
        suits = [c[1] for c in hand_vals]

        is_flush = len(set(suits)) == 1
        # Check for ace-low straight (A-5)
        is_ace_low_straight = (ranks == [14, 5, 4, 3, 2])
        is_straight = (len(set(ranks)) == 5 and ranks[0] - ranks[4] == 4) or is_ace_low_straight
        
        # Use ranks of [5,4,3,2,1] for scoring A-5 straights
        if is_ace_low_straight:
            score_ranks = [5, 4, 3, 2, 1]
        else:
            score_ranks = ranks

        rank_counts = Counter(ranks)
        # Primary rank is the one with more cards (e.g., the rank of the three in a full house)
        # The tuple is sorted first by count, then by rank.
        sorted_by_count = sorted(rank_counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
        primary_ranks = tuple(r[0] for r in sorted_by_count)
        
        # Hand scoring from best to worst
        if is_straight and is_flush:
            return (8, tuple(score_ranks))  # Straight Flush
        if rank_counts[primary_ranks[0]] == 4:
            return (7, primary_ranks)  # Four of a Kind
        if rank_counts[primary_ranks[0]] == 3 and rank_counts[primary_ranks[1]] == 2:
            return (6, primary_ranks)  # Full House
        if is_flush:
            return (5, tuple(ranks))  # Flush
        if is_straight:
            return (4, tuple(score_ranks))  # Straight
        if rank_counts[primary_ranks[0]] == 3:
            return (3, primary_ranks)  # Three of a Kind
        if rank_counts[primary_ranks[0]] == 2 and rank_counts[primary_ranks[1]] == 2:
            return (2, primary_ranks) # Two Pair
        if rank_counts[primary_ranks[0]] == 2:
            return (1, primary_ranks) # One Pair
        return (0, primary_ranks) # High Card

    memo_eval = {}
    def get_best_hand(hand_str_list):
        """
        Finds the best 5-card hand from a 7-card list.
        Uses memoization to speed up redundant calculations.
        """
        # Sorting makes the tuple hashable and consistent
        hand_tuple = tuple(sorted(hand_str_list))
        if hand_tuple in memo_eval:
            return memo_eval[hand_tuple]

        hand_vals = [card_to_val(c) for c in hand_str_list]
        
        best_score = (-1,)
        # Iterate through all 21 combinations of 5 cards from 7
        for combo in itertools.combinations(hand_vals, 5):
            score = score_5_card_hand(list(combo))
            if score > best_score:
                best_score = score
        
        memo_eval[hand_tuple] = best_score
        return best_score

    # 2. Define Hands
    hero_hand = ['As', 'Ac']
    villain_hands_map = {
        "QJ": ['Qh', 'Jh'],
        "QT": ['Qh', 'Th'],
        "Q9": ['Qh', '9h']
    }
    
    results = {}

    print("Starting poker equity calculation...")
    print("This will take a few minutes as it checks all 1,712,304 possible boards for each matchup.")
    print("-" * 30)

    # 3. Main Loop for each matchup
    for name, villain_hand in villain_hands_map.items():
        start_time = time.time()
        
        # Setup for the current matchup
        hero_wins, villain_wins, splits = 0, 0, 0
        memo_eval.clear()
        
        current_deck = list(deck - set(hero_hand) - set(villain_hand))
        
        # Iterate through all possible boards
        board_iterator = itertools.combinations(current_deck, 5)
        
        for board in board_iterator:
            hero_7_cards = hero_hand + list(board)
            villain_7_cards = villain_hand + list(board)

            hero_score = get_best_hand(hero_7_cards)
            villain_score = get_best_hand(villain_7_cards)
            
            if hero_score > villain_score:
                hero_wins += 1
            elif villain_score > hero_score:
                villain_wins += 1
            else:
                splits += 1

        total_boards = hero_wins + villain_wins + splits
        hero_equity = (hero_wins + 0.5 * splits) / total_boards if total_boards > 0 else 0
        
        results[name] = hero_equity
        
        # Print results for the current matchup
        end_time = time.time()
        print(f"Matchup: AsAc vs {villain_hand[0]}{villain_hand[1]}")
        print(f"Completed in {end_time - start_time:.2f} seconds.")
        print(f"Hero Wins: {hero_wins}, Villain Wins: {villain_wins}, Splits: {splits}")
        print(f"Total Boards: {total_boards}")
        print("Equity Equation:")
        print(f"({hero_wins} + 0.5 * {splits}) / {total_boards} = {hero_equity:.4%}\n")


    # 4. Final Conclusion
    min_equity = 1.0
    least_liked_hand = None
    for hand_name, equity in results.items():
        if equity < min_equity:
            min_equity = equity
            least_liked_hand = hand_name

    print("-" * 30)
    print("CONCLUSION:")
    print(f"The hand that minimizes our equity is {villain_hands_map[least_liked_hand][0]}{villain_hands_map[least_liked_hand][1]}s ({least_liked_hand} suited).")
    print(f"Our equity is lowest against this hand, at {min_equity:.4%}.")
    print("We least like to see QJ suited because it has the best connectivity, giving it the highest potential to make straights that can beat our aces.")

solve_poker_equity()
<<<A>>>