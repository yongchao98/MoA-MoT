import collections
from itertools import combinations

def get_hand_rank_name(rank_tuple):
    """Returns the string name of a hand rank tuple."""
    rank = rank_tuple[0]
    if rank == 9: return "Straight Flush"
    if rank == 8: return "Four of a Kind"
    if rank == 7: return "Full House"
    if rank == 6: return "Flush"
    if rank == 5: return "Straight"
    if rank == 4: return "Three of a Kind"
    if rank == 3: return "Two Pair"
    if rank == 2: return "One Pair"
    return "High Card"

def evaluate_hand(cards):
    """
    Evaluates a 7-card hand and returns a tuple representing its rank.
    The tuple is comparable, so higher is better.
    Rank mapping:
    9: Straight Flush
    8: Four of a Kind
    7: Full House
    6: Flush
    5: Straight
    4: Three of a Kind
    3: Two Pair
    2: One Pair
    1: High Card
    """
    if not isinstance(cards, list):
        cards = list(cards)
        
    ranks = sorted([card[0] for card in cards], reverse=True)
    suits = [card[1] for card in cards]
    
    # Check for Flush and Straight Flush
    suit_counts = collections.Counter(suits)
    flush_suit = None
    for suit, count in suit_counts.items():
        if count >= 5:
            flush_suit = suit
            break
            
    if flush_suit:
        flush_cards_ranks = sorted([c[0] for c in cards if c[1] == flush_suit], reverse=True)
        # Check for straight in the flush cards
        # Add Ace as 1 for A-5 straights
        low_straight_ranks = list(flush_cards_ranks)
        if 14 in low_straight_ranks:
            low_straight_ranks.append(1)
        
        # Check for 5 consecutive ranks
        unique_flush_ranks = sorted(list(set(low_straight_ranks)), reverse=True)
        for i in range(len(unique_flush_ranks) - 4):
            if unique_flush_ranks[i] - 4 == unique_flush_ranks[i+4]:
                return (9, unique_flush_ranks[i]) # Straight Flush
        
        return (6, tuple(flush_cards_ranks[:5])) # Flush

    # Check for Straight
    # Add Ace as 1 for A-5 straights
    low_straight_ranks = list(ranks)
    if 14 in low_straight_ranks:
        low_straight_ranks.append(1)
        
    unique_ranks = sorted(list(set(low_straight_ranks)), reverse=True)
    for i in range(len(unique_ranks) - 4):
        if unique_ranks[i] - 4 == unique_ranks[i+4]:
            return (5, unique_ranks[i]) # Straight

    # Check for Quads, Full House, Trips, Pairs
    rank_counts = collections.Counter(ranks)
    # Sort by count desc, then rank desc
    sorted_rank_counts = sorted(rank_counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
    
    counts = [item[1] for item in sorted_rank_counts]
    vals = [item[0] for item in sorted_rank_counts]
    
    if counts[0] == 4: # Four of a Kind
        return (8, vals[0], vals[1])
    if counts[0] == 3 and counts[1] >= 2: # Full House
        return (7, vals[0], vals[1])
    if counts[0] == 3: # Three of a Kind
        kickers = [v for v in vals[1:]][:2]
        return (4, vals[0], tuple(kickers))
    if counts[0] == 2 and counts[1] == 2: # Two Pair
        kicker = vals[2]
        return (3, vals[0], vals[1], kicker)
    if counts[0] == 2: # One Pair
        kickers = [v for v in vals[1:]][:3]
        return (2, vals[0], tuple(kickers))
    
    # High Card
    return (1, tuple(vals[:5]))


def run_simulation(hero_hand, villain_hand):
    """
    Runs a full enumeration of all possible boards for a given matchup.
    """
    # Card values: 2-10, J=11, Q=12, K=13, A=14
    # Suits: s, c, h, d
    deck = []
    for rank in range(2, 15):
        for suit in ['s', 'c', 'h', 'd']:
            deck.append((rank, suit))

    # Remove known cards from deck
    preflop_cards = hero_hand + villain_hand
    remaining_deck = [card for card in deck if card not in preflop_cards]
    
    hero_wins = 0
    villain_wins = 0
    ties = 0

    all_boards = combinations(remaining_deck, 5)
    total_boards = 0

    for board in all_boards:
        total_boards += 1
        hero_full_hand = hero_hand + list(board)
        villain_full_hand = villain_hand + list(board)
        
        hero_rank = evaluate_hand(hero_full_hand)
        villain_rank = evaluate_hand(villain_full_hand)
        
        if hero_rank > villain_rank:
            hero_wins += 1
        elif villain_rank > hero_rank:
            villain_wins += 1
        else:
            ties += 1
            
    hero_equity = (hero_wins + ties / 2) / total_boards
    
    hero_hand_str = "AsAc"
    villain_hand_str_map = {11: "J", 10: "T", 9: "9"}
    villain_hand_str = f"Qd{villain_hand_str_map[villain_hand[1][0]]}d"
    
    print(f"Matchup: {hero_hand_str} vs {villain_hand_str}")
    print(f"Total possible boards: {total_boards}")
    print(f"Aces Wins: {hero_wins}, Opponent Wins: {villain_wins}, Ties: {ties}")
    print(f"Equity Calculation: ({hero_wins} + {ties} / 2) / {total_boards}")
    print(f"Aces Equity: {hero_equity:.4%}\n")
    
    return hero_equity

if __name__ == '__main__':
    # Define hands
    # Ranks: J=11, Q=12, K=13, A=14
    # Suits: s=spades, c=clubs, h=hearts, d=diamonds
    aces = [(14, 's'), (14, 'c')]
    qj_suited = [(12, 'd'), (11, 'd')]
    qt_suited = [(12, 'd'), (10, 'd')]
    q9_suited = [(12, 'd'), (9, 'd')]
    
    opponents = {
        "QJ": qj_suited,
        "QT": qt_suited,
        "Q9": q9_suited
    }
    
    equities = {}
    
    for name, hand in opponents.items():
        equities[name] = run_simulation(aces, hand)
        
    min_equity_hand = min(equities, key=equities.get)
    
    print("-" * 30)
    print(f"The hand that minimizes the aces' equity is {min_equity_hand} suited.")
    print("This is the hand you would least like to see.")
