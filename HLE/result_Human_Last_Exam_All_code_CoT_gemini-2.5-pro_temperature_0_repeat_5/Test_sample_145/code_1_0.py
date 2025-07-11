import random
import itertools
from collections import Counter

def get_card_name(card):
    """Returns the string representation of a card."""
    ranks = {14: 'A', 13: 'K', 12: 'Q', 11: 'J', 10: 'T', 9: '9', 8: '8', 7: '7', 6: '6', 5: '5', 4: '4', 3: '3', 2: '2'}
    suits = {'s': '♠', 'c': '♣', 'h': '♥', 'd': '♦'}
    return f"{ranks[card[0]]}{suits[card[1]]}"

def evaluate_hand(seven_cards):
    """
    Evaluates a 7-card hand and returns the best 5-card hand's score.
    The score is a tuple, allowing direct comparison. Higher tuples are better hands.
    """
    best_score = (-1,)

    for hand in itertools.combinations(seven_cards, 5):
        ranks = sorted([c[0] for c in hand], reverse=True)
        suits = {c[1] for c in hand}
        
        # Handle Ace-low straight (A, 5, 4, 3, 2)
        is_straight = (ranks[0] - ranks[4] == 4 and len(set(ranks)) == 5) or (ranks == [14, 5, 4, 3, 2])
        is_flush = len(suits) == 1
        
        score_ranks = ranks
        if ranks == [14, 5, 4, 3, 2]: # Ace-low straight ranks for scoring
            score_ranks = [5, 4, 3, 2, 1]

        # Score calculation
        score = (-1,)
        counts = Counter(ranks)
        sorted_counts = sorted(counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
        
        # Hand ranking from Straight Flush to High Card
        if is_straight and is_flush:
            score = (8, score_ranks[0]) # Straight Flush
        elif sorted_counts[0][1] == 4:
            score = (7, sorted_counts[0][0], sorted_counts[1][0]) # Four of a Kind
        elif sorted_counts[0][1] == 3 and sorted_counts[1][1] == 2:
            score = (6, sorted_counts[0][0], sorted_counts[1][0]) # Full House
        elif is_flush:
            score = (5, tuple(ranks)) # Flush
        elif is_straight:
            score = (4, score_ranks[0]) # Straight
        elif sorted_counts[0][1] == 3:
            kickers = tuple(r for r in ranks if r != sorted_counts[0][0])
            score = (3, sorted_counts[0][0], kickers) # Three of a Kind
        elif sorted_counts[0][1] == 2 and sorted_counts[1][1] == 2:
            kicker = tuple(r for r in ranks if r != sorted_counts[0][0] and r != sorted_counts[1][0])
            score = (2, sorted_counts[0][0], sorted_counts[1][0], kicker[0]) # Two Pair
        elif sorted_counts[0][1] == 2:
            kickers = tuple(r for r in ranks if r != sorted_counts[0][0])
            score = (1, sorted_counts[0][0], kickers) # One Pair
        else:
            score = (0, tuple(ranks)) # High Card

        if score > best_score:
            best_score = score
            
    return best_score

def run_simulation(hero_hand, villain_hand, num_simulations):
    """Runs a Monte Carlo simulation for a given matchup."""
    ranks = range(2, 15)
    suits = ['s', 'c', 'h', 'd']
    full_deck = list(itertools.product(ranks, suits))

    # Remove known cards from the deck
    deck = [card for card in full_deck if card not in hero_hand and card not in villain_hand]
    
    hero_wins = 0
    villain_wins = 0
    ties = 0

    for _ in range(num_simulations):
        random.shuffle(deck)
        board = deck[:5]
        
        hero_score = evaluate_hand(hero_hand + board)
        villain_score = evaluate_hand(villain_hand + board)

        if hero_score > villain_score:
            hero_wins += 1
        elif villain_score > hero_score:
            villain_wins += 1
        else:
            ties += 1
            
    hero_equity = (hero_wins + ties / 2) / num_simulations
    return hero_equity

# --- Main Execution ---
if __name__ == "__main__":
    # Define hands
    # Ranks: 14=A, 13=K, 12=Q, 11=J, 10=T
    # Suits: s=spades, c=clubs, h=hearts, d=diamonds
    hero_hand = [(14, 's'), (14, 'c')] # Black Aces (AsAc)
    
    villain_hands = {
        "QJ": [(12, 'h'), (11, 'h')], # Red QJ suited (QhJh)
        "QT": [(12, 'h'), (10, 'h')], # Red QT suited (QhTh)
        "Q9": [(12, 'h'), (9, 'h')],  # Red Q9 suited (Qh9h)
    }
    
    num_simulations = 50000 # Increase for more accuracy, decrease for speed
    results = {}

    print(f"Simulating {get_card_name(hero_hand[0])}{get_card_name(hero_hand[1])} vs various hands ({num_simulations} iterations each)...\n")

    for name, hand in villain_hands.items():
        equity = run_simulation(hero_hand, hand, num_simulations)
        results[name] = equity
        villain_hand_str = f"{get_card_name(hand[0])}{get_card_name(hand[1])}"
        print(f"Equity for {get_card_name(hero_hand[0])}{get_card_name(hero_hand[1])} vs. {villain_hand_str}: {equity:.2%}")

    # Find the hand that minimizes the hero's equity
    worst_opponent = min(results, key=results.get)
    min_equity = results[worst_opponent]

    print("\n--- Conclusion ---")
    print(f"The hand that minimizes your equity is {worst_opponent} suited.")
    print(f"Against {worst_opponent} suited, your aces have an equity of approximately {min_equity:.2%}, which is the lowest of the three options.")
    print("This is because more connected cards (like QJ) have a higher probability of making straights and straight flushes, which are the premium hands needed to beat a pair of aces.")
    
    # Determine the final answer choice
    answer_map = {"QJ": "A", "QT": "B", "Q9": "C"}
    final_answer = answer_map.get(worst_opponent, "Error")
    
    print(f"\nFinal Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")
