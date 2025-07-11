import itertools
import random
from collections import Counter

def solve_poker_equity():
    """
    This function simulates the equity of Aces against three different hands
    to determine which one is the most threatening.
    """

    # A robust hand evaluator to determine the rank of a 5-card poker hand.
    # Returns a tuple where a higher value means a better hand.
    def evaluate_5_card_hand(hand):
        ranks = sorted([card[0] for card in hand], reverse=True)
        suits = {card[1] for card in hand}
        is_flush = len(suits) == 1
        is_straight = (len(set(ranks)) == 5 and ranks[0] - ranks[4] == 4) or (ranks == [14, 5, 4, 3, 2]) # Ace-low straight

        # Score is a tuple: (hand_rank, tiebreaker_ranks...)
        if is_straight and is_flush:
            # For wheel straight flush (A-5), the high card is 5.
            score = (8, ranks[0] if ranks != [14, 5, 4, 3, 2] else 5)
        elif 4 in Counter(ranks).values():
            quad_rank = [r for r, c in Counter(ranks).items() if c == 4][0]
            kicker = [r for r, c in Counter(ranks).items() if c == 1][0]
            score = (7, quad_rank, kicker)
        elif sorted(Counter(ranks).values()) == [2, 3]:
            trip_rank = [r for r, c in Counter(ranks).items() if c == 3][0]
            pair_rank = [r for r, c in Counter(ranks).items() if c == 2][0]
            score = (6, trip_rank, pair_rank)
        elif is_flush:
            score = (5, tuple(ranks))
        elif is_straight:
            score = (4, ranks[0] if ranks != [14, 5, 4, 3, 2] else 5)
        elif 3 in Counter(ranks).values():
            trip_rank = [r for r, c in Counter(ranks).items() if c == 3][0]
            kickers = sorted([r for r, c in Counter(ranks).items() if c == 1], reverse=True)
            score = (3, trip_rank, kickers[0], kickers[1])
        elif list(Counter(ranks).values()).count(2) == 2:
            pair_ranks = sorted([r for r, c in Counter(ranks).items() if c == 2], reverse=True)
            kicker = [r for r, c in Counter(ranks).items() if c == 1][0]
            score = (2, pair_ranks[0], pair_ranks[1], kicker)
        elif 2 in Counter(ranks).values():
            pair_rank = [r for r, c in Counter(ranks).items() if c == 2][0]
            kickers = sorted([r for r, c in Counter(ranks).items() if c == 1], reverse=True)
            score = (1, pair_rank, kickers[0], kickers[1], kickers[2])
        else:
            score = (0, tuple(ranks))
        return score

    # Given 7 cards, find the best 5-card hand.
    def get_best_hand(seven_cards):
        return max(itertools.combinations(seven_cards, 5), key=evaluate_5_card_hand)

    # Given best 5-card hands, determine the winner.
    def get_winner(hero_best_5, villain_best_5):
        hero_score = evaluate_5_card_hand(hero_best_5)
        villain_score = evaluate_5_card_hand(villain_best_5)
        if hero_score > villain_score:
            return "hero"
        elif villain_score > hero_score:
            return "villain"
        else:
            return "tie"

    # --- Simulation Setup ---
    ranks_map = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}
    deck = list(itertools.product(ranks_map.values(), ['h', 'd', 'c', 's']))
    
    hero_hand = [(ranks_map['A'], 's'), (ranks_map['A'], 'c')] # Black Aces
    villain_hands_to_test = {
        "QJ": [(ranks_map['Q'], 'h'), (ranks_map['J'], 'h')],
        "QT": [(ranks_map['Q'], 'h'), (ranks_map['T'], 'h')],
        "Q9": [(ranks_map['Q'], 'h'), (ranks_map['9'], 'h')],
    }
    
    num_simulations = 75000
    results = {}
    
    print("Plan: Simulate a large number of hands to calculate the equity (win probability) for Aces against each opponent.")
    print("The hand that gives the Aces the lowest equity is the one you 'least like to see'.\n")

    for name, villain_hand in villain_hands_to_test.items():
        # Setup for this specific matchup
        remaining_deck = [card for card in deck if card not in hero_hand and card not in villain_hand]
        wins = 0
        ties = 0

        # Run Monte Carlo simulation
        for _ in range(num_simulations):
            random.shuffle(remaining_deck)
            board = remaining_deck[:5]
            
            hero_cards_7 = hero_hand + board
            villain_cards_7 = villain_hand + board

            hero_best_hand = get_best_hand(hero_cards_7)
            villain_best_hand = get_best_hand(villain_cards_7)
            
            winner = get_winner(hero_best_hand, villain_best_hand)
            if winner == "hero":
                wins += 1
            elif winner == "tie":
                ties += 1

        equity = (wins + ties / 2) / num_simulations
        results[name] = equity

    print("--- Equity Calculation Results ---")
    min_equity = 1.0
    min_equity_hand = ""
    for hand_name, equity in results.items():
        print(f"Equity for Aces vs. suited {hand_name}: {equity:.2%}")
        if equity < min_equity:
            min_equity = equity
            min_equity_hand = hand_name
    
    print("\n--- Conclusion ---")
    print(f"The hand that minimizes the equity for the Aces is suited {min_equity_hand}.")
    print("This is because the more connected a hand is (like QJ), the greater its potential to make a straight, which makes it more dangerous to an overpair like Aces.")

solve_poker_equity()
<<<A>>>