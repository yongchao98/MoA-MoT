import random
import itertools

def get_hand_rank(hand):
    """
    Evaluates a 5-card hand and returns a numerical rank.
    Higher rank is a better hand. Rank is a tuple for tie-breaking.
    """
    ranks = sorted([card[0] for card in hand], reverse=True)
    suits = [card[1] for card in hand]
    is_flush = len(set(suits)) == 1
    
    # Check for straight (and handle Ace-low straight A-2-3-4-5)
    is_straight = False
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    # Standard straight check
    if len(unique_ranks) >= 5:
        for i in range(len(unique_ranks) - 4):
            if unique_ranks[i] - unique_ranks[i+4] == 4:
                is_straight = True
                straight_high_card = unique_ranks[i]
                break
    # Ace-low straight check
    if not is_straight and set(ranks) == {14, 5, 4, 3, 2}:
        is_straight = True
        straight_high_card = 5
    
    # Get rank counts for pairs, trips, etc.
    rank_counts = {rank: ranks.count(rank) for rank in ranks}
    counts = sorted(rank_counts.values(), reverse=True)
    
    # Determine hand type and return rank tuple
    if is_straight and is_flush: return (8, straight_high_card)
    if counts[0] == 4:
        quad_rank = [r for r, c in rank_counts.items() if c == 4][0]
        kicker = [r for r, c in rank_counts.items() if c == 1][0]
        return (7, quad_rank, kicker)
    if counts == [3, 2]:
        trip_rank = [r for r, c in rank_counts.items() if c == 3][0]
        pair_rank = [r for r, c in rank_counts.items() if c == 2][0]
        return (6, trip_rank, pair_rank)
    if is_flush: return (5, tuple(ranks))
    if is_straight: return (4, straight_high_card)
    if counts[0] == 3:
        trip_rank = [r for r, c in rank_counts.items() if c == 3][0]
        kickers = sorted([r for r in ranks if r != trip_rank], reverse=True)
        return (3, trip_rank, tuple(kickers))
    if counts == [2, 2, 1]:
        pair_ranks = sorted([r for r, c in rank_counts.items() if c == 2], reverse=True)
        kicker = [r for r, c in rank_counts.items() if c == 1][0]
        return (2, tuple(pair_ranks), kicker)
    if counts[0] == 2:
        pair_rank = [r for r, c in rank_counts.items() if c == 2][0]
        kickers = sorted([r for r in ranks if r != pair_rank], reverse=True)
        return (1, pair_rank, tuple(kickers))
    return (0, tuple(ranks))

def get_best_hand(hole_cards, board):
    """Finds the best 5-card hand from a 7-card combination."""
    all_seven_cards = hole_cards + board
    best_rank = (-1,)
    for hand_combination in itertools.combinations(all_seven_cards, 5):
        current_rank = get_hand_rank(list(hand_combination))
        if current_rank > best_rank:
            best_rank = current_rank
    return best_rank

def run_simulation(hero_hand, villain_hand, villain_hand_name, num_simulations=20000):
    """Runs a Monte Carlo simulation for a given matchup."""
    # Ranks: 2-10, J=11, Q=12, K=13, A=14
    # Suits: s=spades, c=clubs, h=hearts, d=diamonds
    deck = list(itertools.product(range(2, 15), ['s', 'h', 'd', 'c']))
    
    # Remove known cards from the deck
    for card in hero_hand + villain_hand:
        deck.remove(card)
        
    hero_wins, ties, villain_wins = 0, 0, 0
    
    for _ in range(num_simulations):
        random.shuffle(deck)
        board = deck[:5]
        
        hero_best = get_best_hand(hero_hand, board)
        villain_best = get_best_hand(villain_hand, board)
        
        if hero_best > villain_best:
            hero_wins += 1
        elif villain_best > hero_best:
            villain_wins += 1
        else:
            ties += 1
            
    total_games = hero_wins + ties + villain_wins
    hero_equity = (hero_wins + ties / 2) / total_games
    
    print(f"Matchup: AsAc vs {villain_hand_name}")
    print(f"  - Results over {total_games} games:")
    print(f"  - Hero Wins: {hero_wins}, Villain Wins: {villain_wins}, Ties: {ties}")
    print(f"  - Hero Equity Equation: ({hero_wins} + {ties} / 2) / {total_games}")
    print(f"  - Hero Equity: {hero_equity:.2%}\n")
    return hero_equity

def main():
    """Main function to run the poker hand analysis."""
    hero_hand = [(14, 's'), (14, 'c')]  # Black Aces (AsAc)

    # Opponent hands (Red Suited)
    villain_hands = {
        "QhJh": [(12, 'h'), (11, 'h')],
        "QhTh": [(12, 'h'), (10, 'h')],
        "Qh9h": [(12, 'h'), (9, 'h')]
    }

    equities = {}
    print("Calculating equity for two black aces against different red suited hands...\n")

    for name, hand in villain_hands.items():
        equities[name] = run_simulation(hero_hand, hand, name)

    min_equity_hand = min(equities, key=equities.get)
    min_equity_value = equities[min_equity_hand]
    
    print("--- Conclusion ---")
    print(f"The hand that minimizes your equity is {min_equity_hand} with an equity of {min_equity_value:.2%}.")
    print("This is because suited connectors have the highest potential to make a straight to beat your aces.")
    print("\nAnswer choice A corresponds to QJ.")

if __name__ == '__main__':
    main()
<<<A>>>