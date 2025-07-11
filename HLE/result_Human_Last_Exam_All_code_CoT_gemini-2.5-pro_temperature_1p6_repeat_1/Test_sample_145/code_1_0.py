import itertools
import random
from collections import Counter

def get_hand_rank(hand):
    """
    Evaluates a 7-card hand and returns a comparable tuple representing its rank.
    Higher ranks correspond to higher tuple values.
    Rank mapping: 9=SF, 8=4K, 7=FH, 6=F, 5=S, 4=3K, 3=2P, 2=1P, 1=HC.
    """
    ranks_map = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}
    
    ranks = sorted([ranks_map[r] for r, s in hand], reverse=True)
    suits = [s for r, s in hand]
    
    suit_counts = Counter(suits)
    flush_suit = next((s for s, c in suit_counts.items() if c >= 5), None)

    flush_ranks = sorted([ranks_map[r] for r, s in hand if s == flush_suit], reverse=True) if flush_suit else []
    
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    
    is_straight = False
    straight_high_card = 0
    # Ace-low straight check (A,2,3,4,5)
    if set([14, 2, 3, 4, 5]).issubset(set(unique_ranks)):
        is_straight = True
        straight_high_card = 5
    for i in range(len(unique_ranks) - 4):
        if unique_ranks[i] - unique_ranks[i+4] == 4:
            is_straight = True
            straight_high_card = max(straight_high_card, unique_ranks[i])
            break
            
    is_straight_flush = False
    sf_high_card = 0
    if flush_suit:
        # Ace-low SF check
        if set([14, 2, 3, 4, 5]).issubset(set(flush_ranks)):
            is_straight_flush = True
            sf_high_card = 5
        # Regular SF check
        for i in range(len(flush_ranks) - 4):
            if flush_ranks[i] - flush_ranks[i+4] == 4:
                is_straight_flush = True
                sf_high_card = max(sf_high_card, flush_ranks[i])
                break

    if is_straight_flush:
        return (9, sf_high_card)

    rank_counts = Counter(ranks)
    counts = sorted(rank_counts.values(), reverse=True)
    ordered_ranks = sorted(rank_counts, key=lambda r: (rank_counts[r], r), reverse=True)

    if counts[0] == 4:
        return (8, ordered_ranks[0], ordered_ranks[1])

    if counts[0] == 3 and counts[1] >= 2:
        return (7, ordered_ranks[0], ordered_ranks[1])

    if flush_suit:
        return (6, tuple(flush_ranks[:5]))

    if is_straight:
        return (5, straight_high_card)

    if counts[0] == 3:
        return (4, ordered_ranks[0], tuple(ordered_ranks[1:3]))
    
    if counts[0] == 2 and counts[1] == 2:
        return (3, ordered_ranks[0], ordered_ranks[1], ordered_ranks[2])

    if counts[0] == 2:
        return (2, ordered_ranks[0], tuple(ordered_ranks[1:4]))

    return (1, tuple(ranks[:5]))

def run_simulation():
    """
    Runs the poker equity simulation and prints the results.
    """
    # Using 's' for spades, 'c' for clubs, 'h' for hearts, 'd' for diamonds
    ranks = '23456789TJQKA'
    suits = 'schd'
    deck = [r + s for r in ranks for s in suits]

    hero_hand = ['As', 'Ac']
    villain_hands_info = {
        'QJ': ['Qh', 'Jh'],
        'QT': ['Qh', 'Th'],
        'Q9': ['Qh', '9h'],
    }
    num_simulations = 200000 # Increase for more accuracy, decrease for speed
    
    results = {}

    print(f"Running {num_simulations} simulations for each matchup...")
    print("-" * 30)

    for name, villain_hand in villain_hands_info.items():
        # Create the remaining deck for simulation
        temp_deck = [card for card in deck if card not in hero_hand and card not in villain_hand]
        
        hero_wins = 0
        villain_wins = 0
        ties = 0

        for _ in range(num_simulations):
            board = random.sample(temp_deck, 5)
            
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
        
        total_runs = hero_wins + villain_wins + ties
        hero_equity = (hero_wins + ties / 2) / total_runs
        results[name] = hero_equity
        
        print(f"Matchup: Black Aces vs. Red {name} suited")
        print(f"Aces wins: {hero_wins}, Villain wins: {villain_wins}, Ties: {ties}")
        print(f"Equity Equation: ({hero_wins} + {ties} / 2) / {total_runs}")
        print(f"Aces Equity: {hero_equity:.2%}")
        print("-" * 30)

    # Find the hand that minimizes the Aces' equity
    worst_opponent = min(results, key=results.get)
    min_equity = results[worst_opponent]

    print("Conclusion:")
    print(f"The hand that minimizes the Aces' equity is {worst_opponent} suited.")
    print(f"Against {worst_opponent} suited, the Aces' equity is only {min_equity:.2%}, the lowest of the three.")
    print("This is because QJ suited has the best connectivity, giving it the highest chance to make a straight against the Aces.")

if __name__ == '__main__':
    run_simulation()
<<<A>>>