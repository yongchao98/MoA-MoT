import itertools
from collections import Counter

def get_hand_strength(seven_cards):
    """
    Evaluates a 7-card hand and returns a comparable tuple representing its strength.
    Hand Rank: 9=StraightFlush, 8=Quads, 7=FullHouse, 6=Flush, 5=Straight,
                 4=Trips, 3=TwoPair, 2=Pair, 1=HighCard.
    The tuple is (hand_rank, kicker1, kicker2, ...).
    """
    ranks = sorted([card[0] for card in seven_cards], reverse=True)
    suits = [card[1] for card in seven_cards]

    # --- Flush and Straight Flush Check ---
    is_flush = False
    flush_suit = None
    suit_counts = Counter(suits)
    for suit, count in suit_counts.items():
        if count >= 5:
            is_flush = True
            flush_suit = suit
            break
    
    flush_ranks = []
    if is_flush:
        flush_ranks = sorted([c[0] for c in seven_cards if c[1] == flush_suit], reverse=True)

    # --- Straight and Straight Flush Check ---
    is_straight = False
    straight_high_rank = -1
    # Use the ranks from a flush if one exists, otherwise all unique ranks
    check_ranks = sorted(list(set(flush_ranks if is_flush else ranks)), reverse=True)
    
    # Ace-low straight check (A,2,3,4,5)
    if set([14, 2, 3, 4, 5]).issubset(check_ranks):
        is_straight = True
        straight_high_rank = 5
    
    for i in range(len(check_ranks) - 4):
        if check_ranks[i] - check_ranks[i+4] == 4:
            is_straight = True
            straight_high_rank = check_ranks[i]
            break

    # --- Return Hand Strength ---
    if is_straight and is_flush:
        return (9, straight_high_rank) # Straight Flush

    # --- Rank-based Hand Evaluation (Quads, Full House, etc.) ---
    rank_counts = Counter(ranks)
    # Sort ranks by count, then by rank value
    sorted_by_count = sorted(rank_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
    
    counts = [item[1] for item in sorted_by_count]
    major_ranks = [item[0] for item in sorted_by_count]

    if counts[0] == 4: # Four of a Kind
        kickers = [r for r in ranks if r != major_ranks[0]]
        return (8, major_ranks[0], kickers[0])
    
    if counts[0] == 3 and counts[1] >= 2: # Full House
        return (7, major_ranks[0], major_ranks[1])

    if is_flush: # Flush
        return (6, tuple(flush_ranks[:5]))

    if is_straight: # Straight
        return (5, straight_high_rank)

    if counts[0] == 3: # Three of a Kind
        kickers = [r for r in ranks if r != major_ranks[0]]
        return (4, major_ranks[0], kickers[0], kickers[1])
    
    if counts[0] == 2 and counts[1] == 2: # Two Pair
        pair_ranks = major_ranks[:2]
        kickers = [r for r in ranks if r not in pair_ranks]
        return (3, pair_ranks[0], pair_ranks[1], kickers[0])

    if counts[0] == 2: # One Pair
        pair_rank = major_ranks[0]
        kickers = [r for r in ranks if r != pair_rank]
        return (2, pair_rank, kickers[0], kickers[1], kickers[2])

    # High Card
    return (1, tuple(ranks[:5]))

def main():
    """Main function to run the equity calculations."""
    # Define ranks and suits for card representation
    # Ranks: T=10, J=11, Q=12, K=13, A=14
    str_ranks = '23456789TJQKA'
    ranks_map = {rank: i for i, rank in enumerate(str_ranks, 2)}
    suits = 'scdh' # spades, clubs, diamonds, hearts

    # --- Setup Hands and Deck ---
    hero_hand_str = ('As', 'Ac') # Black Aces
    villain_hands_str = {
        'QJ suited (QdJd)': ('Qd', 'Jd'),
        'QT suited (QdTd)': ('Qd', 'Td'),
        'Q9 suited (Qd9d)': ('Qd', '9d'),
    }

    # Convert string hands to numerical representation (rank, suit)
    parse = lambda card_str: (ranks_map[card_str[0]], card_str[1])
    hero_hand = tuple(parse(c) for c in hero_hand_str)

    villain_hands = {
        name: tuple(parse(c) for c in hand_tuple)
        for name, hand_tuple in villain_hands_str.items()
    }
    
    full_deck = { (ranks_map[r], s) for r in str_ranks for s in suits }
    
    results = {}
    
    # --- Run Simulation for Each Matchup ---
    for name, villain_hand in villain_hands.items():
        print(f"Calculating equity for {hero_hand_str} vs {name}...")
        
        deck = full_deck - set(hero_hand) - set(villain_hand)
        
        wins, ties, losses = 0, 0, 0
        
        # Iterate through all possible 5-card boards
        board_combinations = itertools.combinations(deck, 5)
        total_boards = 1712304 # C(48, 5)
        
        for board in board_combinations:
            hero_strength = get_hand_strength(hero_hand + board)
            villain_strength = get_hand_strength(villain_hand + board)
            
            if hero_strength > villain_strength:
                wins += 1
            elif villain_strength > hero_strength:
                losses += 1
            else:
                ties += 1

        hero_equity = (wins + ties / 2) / total_boards
        results[name] = hero_equity
        
        print(f"  Hero Wins: {wins}, Villain Wins: {losses}, Ties: {ties}")
        print(f"  Total Boards: {total_boards}")
        print(f"  Equity Equation for Hero: ({wins} + {ties} / 2) / {total_boards}")
        print(f"  Hero Equity: {hero_equity:.3%}\n")
        
    # --- Find the hand we least like to see (lowest equity for us) ---
    least_liked_hand = min(results, key=results.get)
    min_equity = results[least_liked_hand]
    
    print("--- Conclusion ---")
    print("The hand we least like to see is the one that gives our Aces the lowest equity.")
    print(f"The hand that minimizes our equity is: {least_liked_hand}")
    print(f"Against this hand, our equity is {min_equity:.3%}.")


if __name__ == '__main__':
    main()

<<<A>>>