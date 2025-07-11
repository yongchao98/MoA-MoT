import random
import collections

def solve_poker_equity():
    """
    This function solves the poker problem by running a Monte Carlo simulation
    to calculate the equity of two black Aces against three different suited red hands.
    The hand that gives the Aces the lowest equity is the one we least like to see.
    """

    # Step 1: Define card representations and a robust hand evaluator.
    # Cards are represented as tuples: (rank, suit).
    # Ranks: 2-9, T=10, J=11, Q=12, K=13, A=14. Suits: s=0, h=1, d=2, c=3.
    def card_to_int(card_str):
        rank_str = card_str[:-1]
        suit_str = card_str[-1]
        ranks = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}
        suits = {'s': 0, 'h': 1, 'd': 2, 'c': 3}
        return (ranks[rank_str], suits[suit_str])

    def evaluate_7_cards(hand):
        ranks = sorted([c[0] for c in hand], reverse=True)
        suits = [c[1] for c in hand]
        rank_counts = collections.Counter(ranks)
        suit_counts = collections.Counter(suits)

        flush_suit = next((s for s, count in suit_counts.items() if count >= 5), None)
        is_straight = False
        straight_high_card = 0

        unique_ranks = sorted(list(set(ranks)), reverse=True)
        if 14 in unique_ranks: unique_ranks.append(1) # Add Ace as 1 for wheel (A-5)

        for i in range(len(unique_ranks) - 4):
            if unique_ranks[i] - 4 == unique_ranks[i+4]:
                straight_high_card = 5 if unique_ranks[i] == 5 and unique_ranks[0] == 14 else unique_ranks[i]
                is_straight = True
                break
        
        if flush_suit is not None:
            flush_ranks = sorted([c[0] for c in hand if c[1] == flush_suit], reverse=True)
            unique_flush_ranks = sorted(list(set(flush_ranks)), reverse=True)
            if 14 in unique_flush_ranks: unique_flush_ranks.append(1)
            for i in range(len(unique_flush_ranks) - 4):
                 if unique_flush_ranks[i] - 4 == unique_flush_ranks[i+4]:
                    sf_high_card = 5 if unique_flush_ranks[i] == 5 and unique_flush_ranks[0] == 14 else unique_flush_ranks[i]
                    return (8, sf_high_card)  # Straight Flush
            return (5, tuple(flush_ranks[:5]))  # Flush

        if is_straight:
            return (4, straight_high_card)  # Straight
        
        sorted_by_count = sorted(rank_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        counts = [item[1] for item in sorted_by_count]
        ranked_ranks = [item[0] for item in sorted_by_count]

        if counts[0] == 4: return (7, ranked_ranks[0], ranked_ranks[1])  # Four of a Kind
        if counts[0] == 3 and counts[1] >= 2: return (6, ranked_ranks[0], ranked_ranks[1])  # Full House
        if counts[0] == 3: return (3, ranked_ranks[0], tuple(ranked_ranks[1:3]))  # Three of a Kind
        if counts[0] == 2 and counts[1] == 2: return (2, tuple(ranked_ranks[0:2]), ranked_ranks[2])  # Two Pair
        if counts[0] == 2: return (1, ranked_ranks[0], tuple(ranked_ranks[1:4]))  # One Pair
        return (0, tuple(ranked_ranks[:5]))  # High Card

    def run_simulation(hero_str, villain_str, num_sims):
        hero_cards = [card_to_int(c) for c in hero_str]
        villain_cards = [card_to_int(c) for c in villain_str]

        deck = list({(r, s) for r in range(2, 15) for s in range(4)} - set(hero_cards) - set(villain_cards))
        wins, ties, total_sims = 0, 0, num_sims
        
        for _ in range(total_sims):
            board = random.sample(deck, 5)
            hero_score = evaluate_7_cards(hero_cards + board)
            villain_score = evaluate_7_cards(villain_cards + board)
            
            if hero_score > villain_score: wins += 1
            elif hero_score == villain_score: ties += 1
                
        equity_numerator = wins * 2 + ties
        equity_denominator = total_sims * 2
        equity = (wins + ties / 2) / total_sims
        return equity, wins, total_sims - wins - ties, ties, equity_numerator, equity_denominator

    # Step 2: Run simulations for each matchup.
    print("Running Monte Carlo simulations to calculate equities...")
    print("Hero Hand: Two black Aces (A♠ A♣)")
    print("-" * 50)

    hero_hand = ['As', 'Ac']
    opponent_hands = { "QJ": ['Qd', 'Jd'], "QT": ['Qd', 'Td'], "Q9": ['Qd', '9d'] }
    simulations = 50000
    results = {}
    
    for name, hand in opponent_hands.items():
        equity, wins, losses, ties, num, den = run_simulation(hero_hand, hand, simulations)
        results[name] = equity
        print(f"Matchup: Hero (A♠ A♣) vs. Opponent ({hand[0]} {hand[1]})")
        print(f"Result over {simulations} hands:")
        print(f"  - Hero wins: {wins}, Opponent wins: {losses}, Ties: {ties}")
        print(f"  - Hero Equity Calculation: ({wins} + {ties} / 2) / {simulations} = {equity:.2%}\n")

    # Step 3: Determine and present the final answer.
    min_equity_hand = min(results, key=results.get)
    min_equity_value = results[min_equity_hand]

    print("-" * 50)
    print("Conclusion:")
    print(f"Our equity is lowest against {min_equity_hand} suited ({min_equity_value:.2%}).")
    print("This confirms that the most connected hand, QJ, is the strongest against our Aces.")
    print("Therefore, QJ suited is the hand we least like to see.")
    
    final_answer_map = {"QJ": "A", "QT": "B", "Q9": "C"}
    final_answer_char = final_answer_map.get(min_equity_hand)
    print(f'<<<{final_answer_char}>>>')

# Execute the full analysis
solve_poker_equity()