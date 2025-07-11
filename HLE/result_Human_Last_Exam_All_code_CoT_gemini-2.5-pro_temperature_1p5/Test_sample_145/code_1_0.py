import random

def evaluate_hand(hand):
    """
    Evaluates a 7-card poker hand and returns a tuple representing its strength.
    The tuple is structured for easy comparison: (hand_rank, kicker1, kicker2, ...).
    Higher is better.
    Hand Ranks: 8=StraightFlush, 7=Quads, 6=FullHouse, 5=Flush, 4=Straight,
                  3=Trips, 2=TwoPair, 1=Pair, 0=HighCard.
    """
    ranks_str = '23456789TJQKA'
    
    # Parse cards into numeric ranks and suits
    ranks = sorted([ranks_str.find(c[0]) for c in hand], reverse=True)
    suits = [c[1] for c in hand]
    
    # Check for flush
    is_flush = False
    flush_suit = None
    for s in "shdc":
        if suits.count(s) >= 5:
            is_flush = True
            flush_suit = s
            break
            
    # Check for straight
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    is_straight = False
    straight_high_card = -1
    # Check for Ace-low straight (A, 2, 3, 4, 5)
    if set([12, 3, 2, 1, 0]).issubset(set(unique_ranks)):
        is_straight = True
        straight_high_card = 3 # Rank of '5'
    # Check for other straights
    for i in range(len(unique_ranks) - 4):
        if unique_ranks[i] - unique_ranks[i+4] == 4:
            is_straight = True
            straight_high_card = unique_ranks[i]
            break

    # Straight Flush
    if is_flush and is_straight:
        flush_ranks = sorted([ranks_str.find(c[0]) for c in hand if c[1] == flush_suit], reverse=True)
        # Check for ace-low straight flush
        if set([12, 3, 2, 1, 0]).issubset(set(flush_ranks)):
            return (8, 3) # Royal flush is a straight flush to ace. A-5 is to 5
        # Check for other straight flushes
        for i in range(len(flush_ranks) - 4):
            if flush_ranks[i] - flush_ranks[i+4] == 4:
                return (8, flush_ranks[i])

    # Rank counts for Quads, Full House, etc.
    rank_counts = {r: ranks.count(r) for r in unique_ranks}
    # Sort by count desc, then rank desc
    counts = sorted(rank_counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
    
    # Quads (Four of a Kind)
    if counts[0][1] == 4:
        quad_rank = counts[0][0]
        kickers = [r for r in ranks if r != quad_rank]
        return (7, quad_rank, kickers[0])

    # Full House
    if counts[0][1] == 3 and counts[1][1] >= 2:
        return (6, counts[0][0], counts[1][0])

    # Flush
    if is_flush:
        flush_ranks = sorted([ranks_str.find(c[0]) for c in hand if c[1] == flush_suit], reverse=True)
        return (5, tuple(flush_ranks[:5]))
    
    # Straight
    if is_straight:
        return (4, straight_high_card)
        
    # Trips (Three of a Kind)
    if counts[0][1] == 3:
        trip_rank = counts[0][0]
        kickers = [r for r in ranks if r != trip_rank]
        return (3, trip_rank, tuple(kickers[:2]))

    # Two Pair
    if counts[0][1] == 2 and counts[1][1] == 2:
        pair_ranks = sorted([counts[0][0], counts[1][0]], reverse=True)
        kickers = [r for r in ranks if r not in pair_ranks]
        return (2, tuple(pair_ranks), kickers[0])

    # One Pair
    if counts[0][1] == 2:
        pair_rank = counts[0][0]
        kickers = [r for r in ranks if r != pair_rank]
        return (1, pair_rank, tuple(kickers[:3]))
        
    # High Card
    return (0, tuple(ranks[:5]))


def run_simulation(hero_hand, villain_hand, simulations=20000):
    """Runs a Monte Carlo simulation for a given Hold'em matchup."""
    ranks = '23456789TJQKA'
    suits = 'shdc'
    full_deck = [r + s for r in ranks for s in suits]

    # Remove known cards from the deck
    deck = [card for card in full_deck if card not in hero_hand + villain_hand]
    
    hero_wins = 0
    villain_wins = 0
    ties = 0

    for _ in range(simulations):
        random.shuffle(deck)
        board = deck[:5]
        
        hero_full_hand = hero_hand + board
        villain_full_hand = villain_hand + board
        
        hero_score = evaluate_hand(hero_full_hand)
        villain_score = evaluate_hand(villain_full_hand)
        
        if hero_score > villain_score:
            hero_wins += 1
        elif villain_score > hero_score:
            villain_wins += 1
        else:
            ties += 1
            
    hero_equity = (hero_wins + ties / 2) / simulations * 100
    
    print(f"Matchup: {hero_hand} vs {villain_hand}")
    print(f"Simulations: {simulations}")
    print(f"Aces Wins: {hero_wins}, Opponent Wins: {villain_wins}, Ties: {ties}")
    # Show the calculation for equity
    print(f"Aces Equity calculation: ({hero_wins} + {ties}/2) / {simulations} = {hero_equity:.2f}%\n")
    return hero_equity

# --- Main ---
# Our hand: two black aces
hero = ['As', 'Ac']

# Opponent hands: three different red suited hands
villains = {
    "QJ": ['Qh', 'Jh'],
    "QT": ['Qh', 'Th'],
    "Q9": ['Qh', '9h']
}

results = {}
for name, hand in villains.items():
    results[name] = run_simulation(hero, hand)

# Determine which hand is worst for the hero (lowest equity)
worst_hand_for_hero = min(results, key=results.get)
print(f"The hand that minimizes the Aces' equity is {worst_hand_for_hero} suited.")
print("This is the hand you least like to see.")

<<<A>>>