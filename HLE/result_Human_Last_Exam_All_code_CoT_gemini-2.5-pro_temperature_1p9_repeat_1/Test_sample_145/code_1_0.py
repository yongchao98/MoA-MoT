import itertools
from collections import Counter
import random

def get_rank(card):
    """Returns the rank of a card."""
    return card[0]

def evaluate_hand(five_card_hand):
    """
    Evaluates a 5-card hand and returns a comparable score tuple.
    Higher tuple value means a better hand.
    (hand_rank, high_card1, high_card2, ...)
    """
    ranks = sorted([get_rank(c) for c in five_card_hand], reverse=True)
    suits = [c[1] for c in five_card_hand]
    
    is_flush = len(set(suits)) == 1
    is_wheel = (ranks == [14, 5, 4, 3, 2])
    is_straight = all(ranks[i] == ranks[i+1] + 1 for i in range(len(ranks)-1)) or is_wheel

    if is_straight and is_flush:
        return (9, ranks[0] if not is_wheel else 5)  # Straight Flush (A-5 wheel ranked by 5)

    rank_counts = Counter(ranks)
    counts = sorted(rank_counts.values(), reverse=True)
    major_ranks = sorted(rank_counts.keys(), key=lambda r: (rank_counts[r], r), reverse=True)

    if counts[0] == 4:
        return (8, major_ranks[0], major_ranks[1])  # Four of a Kind
    if counts == [3, 2]:
        return (7, major_ranks[0], major_ranks[1])  # Full House
    if is_flush:
        return (6, tuple(ranks))  # Flush
    if is_straight:
        return (5, ranks[0] if not is_wheel else 5)  # Straight
    if counts[0] == 3:
        return (4, major_ranks[0], tuple(major_ranks[1:3]))  # Three of a Kind
    if counts[:2] == [2, 2]:
        return (3, major_ranks[0], major_ranks[1], major_ranks[2])  # Two Pair
    if counts[0] == 2:
        return (2, major_ranks[0], tuple(major_ranks[1:4]))  # One Pair
        
    return (1, tuple(ranks))  # High Card

def find_best_hand(seven_cards):
    """From 7 cards, find the best 5-card hand by checking all 21 combinations."""
    best_score = (0,)
    for hand_combo in itertools.combinations(seven_cards, 5):
        score = evaluate_hand(hand_combo)
        if score > best_score:
            best_score = score
    return score

def calculate_equity_simulation(hero_hand, villain_hand, n_samples=200000):
    """Runs a Monte Carlo simulation to calculate equity."""
    # Ranks: 2-10, J=11, Q=12, K=13, A=14
    # Suits: s=0, c=1 (black), h=2, d=3 (red)
    all_ranks = list(range(2, 15))
    all_suits = list(range(4))
    full_deck = [c for c in itertools.product(all_ranks, all_suits)]

    # Remove dealt cards from the deck
    for card in hero_hand + villain_hand:
        full_deck.remove(card)

    wins, ties, losses = 0, 0, 0
    # Use a fixed seed for reproducible results
    random.seed(42)
    
    for _ in range(n_samples):
        board = random.sample(full_deck, 5)
        
        # find_best_hand requires iterating all 21 combos, which is slow.
        # A fast C-based evaluator or a lookup table is standard.
        # Here we manually find the best hand for both players.
        best_hero_score = max(evaluate_hand(combo) for combo in itertools.combinations(hero_hand + board, 5))
        best_villain_score = max(evaluate_hand(combo) for combo in itertools.combinations(villain_hand + board, 5))

        if best_hero_score > best_villain_score:
            wins += 1
        elif best_hero_score < best_villain_score:
            losses += 1
        else:
            ties += 1

    equity = (wins + ties / 2) / n_samples
    return equity, wins, ties, n_samples

def card_to_str(card):
    """Helper function to print cards in a readable format."""
    rank_str = {10: 'T', 11: 'J', 12: 'Q', 13: 'K', 14: 'A'}.get(card[0], str(card[0]))
    suit_str = {0: 's', 1: 'c', 2: 'h', 3: 'd'}[card[1]]
    return rank_str + suit_str

def main():
    """Main function to run the analysis."""
    # Define hands: Ranks 2-14, Suits 0-3 (s,c,h,d)
    hero_hand = [(14, 0), (14, 1)]  # Black Aces: AsAc
    
    villain_hands = {
        "QJ suited": [(12, 2), (11, 2)],  # Red QJ: QhJh
        "QT suited": [(12, 2), (10, 2)],  # Red QT: QhTh
        "Q9 suited": [(12, 2), (9, 2)],   # Red Q9: Qh9h
    }
    
    hero_str = card_to_str(hero_hand[0]) + card_to_str(hero_hand[1])
    results = {}
    
    print(f"Calculating equity for {hero_str} against different hands.")
    print("This is done by simulating 200,000 random boards for each matchup.\n")

    min_equity = 1.0
    worst_hand_for_aa = None
    
    for name, hand in villain_hands.items():
        villain_str = card_to_str(hand[0]) + card_to_str(hand[1])
        equity, wins, ties, total = calculate_equity_simulation(hero_hand, hand)
        results[name] = equity
        
        if equity < min_equity:
            min_equity = equity
            worst_hand_for_aa = name
            
        print(f"Matchup: {hero_str} vs. {villain_str}")
        print(f"Result over {total} boards:")
        print(f"  - Wins for {hero_str}: {wins}")
        print(f"  - Ties: {ties}")
        print(f"  - Losses for {hero_str}: {total - wins - ties}")
        print(f"Equity Equation: ({wins} + {ties} / 2) / {total}")
        print(f"Final Equity for {hero_str}: {equity:.4f} or {equity:.2%}\n")

    print("--- Conclusion ---")
    print("The hand we least like to see is the one that gives our Aces the lowest equity.")
    print(f"The simulation shows the hand that minimizes {hero_str}'s equity is {worst_hand_for_aa}.")
    print("This is because its superior connectivity gives it a higher chance of making a straight to beat our pair.")

main()
<<<A>>>