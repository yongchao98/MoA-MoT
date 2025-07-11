import itertools
import random
from collections import defaultdict

def get_hand_rank(hand):
    """Calculates the rank of a 5-card hand."""
    if len(hand) != 5:
        raise ValueError("Hand must contain 5 cards.")
    
    values = sorted([int(r) for r in '..23456789TJQKA'.index(c[0])] for c in hand)
    ranks = [v[0] for v in values]
    suits = [c[1] for c in hand]
    
    is_flush = len(set(suits)) == 1
    # Straight check, including A-5
    is_straight = (max(ranks) - min(ranks) == 4 and len(set(ranks)) == 5) or (ranks == [2, 3, 4, 5, 14])

    if is_straight and is_flush:
        return (8, ranks[-1])  # Straight Flush
    
    counts = defaultdict(int)
    for r in ranks:
        counts[r] += 1
    
    # Sort ranks by count, then by rank value
    sorted_ranks = sorted(counts.keys(), key=lambda r: (counts[r], r), reverse=True)
    
    if counts[sorted_ranks[0]] == 4:
        return (7, sorted_ranks[0], sorted_ranks[1])  # Four of a Kind
    
    if counts[sorted_ranks[0]] == 3 and counts[sorted_ranks[1]] == 2:
        return (6, sorted_ranks[0], sorted_ranks[1])  # Full House
        
    if is_flush:
        return (5, tuple(sorted(ranks, reverse=True))) # Flush
        
    if is_straight:
        # Handle A-5 straight high card (5)
        if ranks == [2, 3, 4, 5, 14]:
            return (4, 5)
        return (4, ranks[-1])  # Straight

    if counts[sorted_ranks[0]] == 3:
        kickers = tuple(sorted([r for r in sorted_ranks[1:]], reverse=True))
        return (3, sorted_ranks[0]) + kickers # Three of a Kind

    if counts[sorted_ranks[0]] == 2 and counts[sorted_ranks[1]] == 2:
        kicker = sorted_ranks[2]
        pairs = tuple(sorted([sorted_ranks[0], sorted_ranks[1]], reverse=True))
        return (2, pairs, kicker) # Two Pair

    if counts[sorted_ranks[0]] == 2:
        kickers = tuple(sorted([r for r in sorted_ranks[1:]], reverse=True))
        return (1, sorted_ranks[0]) + kickers # One Pair

    return (0, tuple(sorted(ranks, reverse=True))) # High Card

def get_best_hand(cards):
    """Finds the best 5-card hand from a list of 7 cards."""
    return max(get_hand_rank(combo) for combo in itertools.combinations(cards, 5))

def calculate_equity(hero_hand, villain_hand, num_simulations=10000):
    """Monte Carlo simulation to calculate equity."""
    deck = {r + s for r in '23456789TJQKA' for s in 'shdc'}
    deck -= set(hero_hand)
    deck -= set(villain_hand)
    
    wins = 0
    ties = 0

    for _ in range(num_simulations):
        board = random.sample(list(deck), 5)
        hero_best = get_best_hand(hero_hand + board)
        villain_best = get_best_hand(villain_hand + board)

        if hero_best > villain_best:
            wins += 1
        elif hero_best == villain_best:
            ties += 1
            
    return (wins + ties / 2) / num_simulations

def get_combos(hand_str, hero_hand=[]):
    """Get all card combinations for a hand string like 'AKs' or 'TT'."""
    ranks = '23456789TJQKA'
    suits = 'shdc'
    r1, r2 = hand_str[0], hand_str[1]
    
    # Remove hero's cards from potential combos
    blocked_cards = set(hero_hand)
    
    combos = []
    if r1 == r2: # Pocket Pair
        combos = [tuple(sorted(c)) for c in itertools.combinations([r1 + s for s in suits], 2)]
    elif len(hand_str) == 3 and hand_str[2] == 's': # Suited
        combos = [tuple(sorted((r1 + s, r2 + s))) for s in suits]
    else: # Offsuit
        combos = [tuple(sorted((r1 + s1, r2 + s2))) for s1 in suits for s2 in suits if s1 != s2]

    # Filter out combos that use a blocked card
    valid_combos = [c for c in combos if not (c[0] in blocked_cards or c[1] in blocked_cards)]
    return valid_combos

def main():
    # --- Configuration ---
    hero_options = {
        "QJs": ['Qh', 'Jh'],
        "99": ['9h', '9s'],
        "AJo": ['Ah', 'Jd'],
        "AKo": ['Ah', 'Kd'],
    }
    
    # Villain's tight calling range on the bubble
    villain_range_notation = ['TT', 'JJ', 'QQ', 'KK', 'AA', 'AQs', 'AKo']
    
    results = {}
    print("Calculating showdown equity for each hand vs a tight bubble calling range (TT+, AQs+, AKo)\n")

    for hero_name, hero_cards in hero_options.items():
        print(f"--- Analyzing Hero Hand: {hero_name} ({''.join(hero_cards)}) ---")
        
        total_weighted_equity = 0
        total_combos = 0
        
        equity_details = []

        for hand_notation in villain_range_notation:
            villain_combos = get_combos(hand_notation, hero_hand=hero_cards)
            num_combos = len(villain_combos)
            
            if num_combos == 0:
                continue

            # Use the first valid combo as a representative for the equity calculation
            representative_villain_hand = list(villain_combos[0])
            equity = calculate_equity(hero_cards, representative_villain_hand, num_simulations=2000) # Lower sims for speed
            
            total_weighted_equity += equity * num_combos
            total_combos += num_combos
            equity_details.append({'notation': hand_notation, 'equity': equity, 'combos': num_combos})

        average_equity = total_weighted_equity / total_combos if total_combos > 0 else 0
        results[hero_name] = average_equity
        
        print(f"Average Equity for {hero_name}: {average_equity:.2%}\n")

    # Find and print the best hand
    best_hand = max(results, key=results.get)
    print("--- Conclusion ---")
    print("Based on the simulations, the hand with the highest equity when called is:")
    for hand, equity in sorted(results.items(), key=lambda item: item[1], reverse=True):
        print(f"{hand}: {equity:.2%}")
        
    print(f"\nThe best hand to jam is {best_hand} because it has the highest profitability in a showdown.")

    # Fulfilling request to show the final equation for the best hand
    print(f"\nExample of the final weighted average equity calculation for {best_hand}:")
    hero_cards = hero_options[best_hand]
    equation_parts = []
    total_combos = 0
    total_weighted_equity = 0

    for hand_notation in villain_range_notation:
        villain_combos = get_combos(hand_notation, hero_hand=hero_cards)
        num_combos = len(villain_combos)
        if num_combos > 0:
            equity = calculate_equity(hero_cards, list(villain_combos[0]), num_simulations=2000)
            equation_parts.append(f"({equity:.3f} * {num_combos} combos of {hand_notation})")
            total_combos += num_combos
            total_weighted_equity += equity * num_combos

    final_equity = total_weighted_equity / total_combos
    equation_str = " + ".join(equation_parts)
    print(f"({equation_str}) / {total_combos} total combos = {final_equity:.2%}")


if __name__ == '__main__':
    main()
