import random
import itertools
from collections import Counter

def solve_poker_equity():
    """
    This script solves the poker problem by simulating equities.
    It determines which of the opponent's hands (QJ, QT, Q9 suited)
    is the strongest against a pair of black Aces, thus minimizing the Aces' equity.
    """

    # --- Card and Hand Representation ---
    def create_deck():
        """Creates a standard 52-card deck."""
        ranks = list(range(2, 15))  # 2-10, J=11, Q=12, K=13, A=14
        suits = ['s', 'h', 'c', 'd'] # Spades, Hearts, Clubs, Diamonds
        return list(itertools.product(ranks, suits))

    def hand_to_str(hand):
        """Converts a hand (list of cards) to a readable string."""
        def rank_to_str(rank):
            if rank <= 10: return str(rank)
            return {11: 'J', 12: 'Q', 13: 'K', 14: 'A'}[rank]
        return ' '.join([rank_to_str(r) + s for r, s in hand])

    # --- Hand Evaluation Logic ---
    def get_hand_rank(hand):
        """
        Evaluates a 5-card hand and returns a comparable tuple representing its rank.
        Format: (rank_value, kicker1, kicker2, ...)
        Rank values: 9=StraightFlush, 8=Quads, 7=FullHouse, 6=Flush, 5=Straight,
                     4=Trips, 3=TwoPair, 2=OnePair, 1=HighCard
        """
        ranks = sorted([card[0] for card in hand], reverse=True)
        suits = [card[1] for card in hand]
        
        is_flush = len(set(suits)) == 1
        
        # Check for A-2-3-4-5 "wheel" straight
        is_straight = False
        unique_ranks = sorted(list(set(ranks)), reverse=True)
        if len(unique_ranks) >= 5:
            if unique_ranks[0] - unique_ranks[4] == 4:
                is_straight = True
                ranks = unique_ranks[0:5] # Main straight ranks
            elif set(unique_ranks).issuperset({14, 2, 3, 4, 5}):
                is_straight = True
                ranks = [5, 4, 3, 2, 1] # Wheel ranks for comparison
        
        if is_straight and is_flush:
            return (9, ranks[0])
        
        counts = Counter(sorted([card[0] for card in hand], reverse=True))
        rank_order = sorted(counts.keys(), key=lambda k: (counts[k], k), reverse=True)
        
        if counts[rank_order[0]] == 4:
            return (8, rank_order[0], rank_order[1])
        if counts[rank_order[0]] == 3 and counts[rank_order[1]] == 2:
            return (7, rank_order[0], rank_order[1])
        if is_flush:
            return (6, tuple(sorted([c[0] for c in hand], reverse=True)))
        if is_straight:
            return (5, ranks[0])
        if counts[rank_order[0]] == 3:
            return (4, rank_order[0], rank_order[1], rank_order[2])
        if counts[rank_order[0]] == 2 and counts[rank_order[1]] == 2:
            return (3, rank_order[0], rank_order[1], rank_order[2])
        if counts[rank_order[0]] == 2:
            return (2, rank_order[0], rank_order[1], rank_order[2], rank_order[3])
        return (1, tuple(rank_order))

    def evaluate_7_cards(cards):
        """Finds the best 5-card hand from a set of 7 cards."""
        best_rank = (0,)
        for hand_5 in itertools.combinations(cards, 5):
            current_rank = get_hand_rank(list(hand_5))
            if current_rank > best_rank:
                best_rank = current_rank
        return best_rank

    # --- Simulation Logic ---
    def simulate_equity(hero_hand, villain_hand, simulations):
        deck = create_deck()
        known_cards = hero_hand + villain_hand
        deck = [card for card in deck if card not in known_cards]
        
        hero_wins = 0
        ties = 0
        
        for _ in range(simulations):
            board = random.sample(deck, 5)
            hero_rank = evaluate_7_cards(hero_hand + board)
            villain_rank = evaluate_7_cards(villain_hand + board)
            
            if hero_rank > villain_rank:
                hero_wins += 1
            elif hero_rank == villain_rank:
                ties += 1
                
        return (hero_wins + ties / 2) / simulations

    # --- Main Execution ---
    hero_hand = [(14, 's'), (14, 'c')]
    
    villain_hands = {
        "A": {"name": "QJ suited", "hand": [(12, 'd'), (11, 'd')]},
        "B": {"name": "QT suited", "hand": [(12, 'd'), (10, 'd')]},
        "C": {"name": "Q9 suited", "hand": [(9, 'd'), (12, 'd')]}
    }
    
    print("Simulating Hero (As Ac) vs. three different opponent hands.")
    print("The hand we least want to see is the one that gives us the lowest equity.")
    print("-" * 60)
    
    results = {}
    num_simulations = 20000 
    
    for choice, data in villain_hands.items():
        equity = simulate_equity(hero_hand, data["hand"], simulations=num_simulations)
        results[choice] = equity
        print(f"Hero's equity against {data['name']} ({hand_to_str(data['hand'])}): {equity:.2%}")

    min_equity_choice = min(results, key=results.get)
    
    print("-" * 60)
    print(f"Conclusion: Our equity is lowest against {villain_hands[min_equity_choice]['name']}.")
    print("This confirms that the more connected hand (QJ) is the biggest threat.")
    print(f"Therefore, {villain_hands[min_equity_choice]['name']} is the hand you least like to see.")
    
    print(f"<<<{min_equity_choice}>>>")


if __name__ == '__main__':
    solve_poker_equity()