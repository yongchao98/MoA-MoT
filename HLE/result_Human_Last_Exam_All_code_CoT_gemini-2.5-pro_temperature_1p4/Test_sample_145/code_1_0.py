import itertools
from collections import Counter

def create_deck():
    """Creates a standard 52-card deck."""
    ranks = list(range(2, 15))  # 2-10, J(11), Q(12), K(13), A(14)
    suits = ['s', 'c', 'h', 'd'] # Spades, Clubs, Hearts, Diamonds
    return list(itertools.product(ranks, suits))

def parse_card(card_str):
    """Converts a card string like 'As' to a tuple (14, 's')."""
    rank_map = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8,
                '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}
    rank = rank_map[card_str[0]]
    suit = card_str[1]
    return (rank, suit)

def get_hand_rank(seven_cards):
    """
    Evaluates the best 5-card hand from a set of 7 cards.
    Returns a tuple for comparison, e.g., (8, 14) for a royal flush.
    Hand ranks: 8=StraightFlush, 7=Quads, 6=FullHouse, 5=Flush, 4=Straight,
                3=Trips, 2=TwoPair, 1=Pair, 0=HighCard.
    """
    ranks = sorted([card[0] for card in seven_cards], reverse=True)
    suits = [card[1] for card in seven_cards]
    
    # Check for Flush
    suit_counts = Counter(suits)
    flush_suit = None
    if suit_counts.most_common(1)[0][1] >= 5:
        flush_suit = suit_counts.most_common(1)[0][0]

    # Check for Straight
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    
    def find_straight(rank_list):
        if 14 in rank_list: # Check for A-2-3-4-5 straight
            rank_list_with_low_ace = rank_list + [1]
        else:
            rank_list_with_low_ace = rank_list
            
        if len(rank_list_with_low_ace) >= 5:
            # Check for 5 consecutive cards
            for i in range(len(rank_list_with_low_ace) - 4):
                is_consecutive = True
                for j in range(4):
                    if rank_list_with_low_ace[i+j] - 1 != rank_list_with_low_ace[i+j+1]:
                        is_consecutive = False
                        break
                if is_consecutive:
                    return rank_list_with_low_ace[i] # Return high card of straight
        return None

    # Evaluate hand strength from highest to lowest
    
    # 1. Straight Flush
    if flush_suit:
        flush_ranks = sorted([c[0] for c in seven_cards if c[1] == flush_suit], reverse=True)
        straight_flush_high = find_straight(flush_ranks)
        if straight_flush_high is not None:
            return (8, straight_flush_high)

    # 2. Four of a Kind
    rank_counts = Counter(ranks)
    if rank_counts.most_common(1)[0][1] == 4:
        quad_rank = rank_counts.most_common(1)[0][0]
        kicker = max(r for r in ranks if r != quad_rank)
        return (7, quad_rank, kicker)

    # 3. Full House
    if rank_counts.most_common(2)[0][1] == 3 and rank_counts.most_common(2)[1][1] >= 2:
        trips_rank = rank_counts.most_common(1)[0][0]
        # Find the highest pair among remaining cards
        pair_ranks = [r for r,c in rank_counts.items() if c >= 2 and r != trips_rank]
        pair_rank = max(pair_ranks)
        return (6, trips_rank, pair_rank)

    # 4. Flush
    if flush_suit:
        flush_ranks = sorted([c[0] for c in seven_cards if c[1] == flush_suit], reverse=True)
        return (5, tuple(flush_ranks[:5]))
    
    # 5. Straight
    straight_high = find_straight(unique_ranks)
    if straight_high is not None:
        return (4, straight_high)

    # 6. Three of a Kind
    if rank_counts.most_common(1)[0][1] == 3:
        trips_rank = rank_counts.most_common(1)[0][0]
        kickers = [r for r in ranks if r != trips_rank]
        return (3, trips_rank, tuple(kickers[:2]))

    # 7. Two Pair
    if rank_counts.most_common(2)[0][1] == 2 and rank_counts.most_common(2)[1][1] == 2:
        pairs = rank_counts.most_common(2)
        high_pair = max(pairs[0][0], pairs[1][0])
        low_pair = min(pairs[0][0], pairs[1][0])
        kicker = max(r for r in ranks if r != high_pair and r != low_pair)
        return (2, high_pair, low_pair, kicker)

    # 8. One Pair
    if rank_counts.most_common(1)[0][1] == 2:
        pair_rank = rank_counts.most_common(1)[0][0]
        kickers = [r for r in ranks if r != pair_rank]
        return (1, pair_rank, tuple(kickers[:3]))
    
    # 9. High Card
    return (0, tuple(ranks[:5]))

def calculate_equity(hero_hand, villain_hand):
    """
    Calculates hero's equity against a villain by iterating all possible boards.
    """
    deck = create_deck()
    
    # Remove known cards from deck
    deck.remove(hero_hand[0])
    deck.remove(hero_hand[1])
    deck.remove(villain_hand[0])
    deck.remove(villain_hand[1])

    wins = 0
    ties = 0
    total_boards = 0
    
    for board in itertools.combinations(deck, 5):
        total_boards += 1
        hero_7_cards = hero_hand + list(board)
        villain_7_cards = villain_hand + list(board)
        
        hero_rank = get_hand_rank(hero_7_cards)
        villain_rank = get_hand_rank(villain_7_cards)
        
        if hero_rank > villain_rank:
            wins += 1
        elif hero_rank == villain_rank:
            ties += 1
            
    equity = (wins + ties / 2) / total_boards
    return equity, wins, ties, total_boards

if __name__ == "__main__":
    hero_hand_str = ('As', 'Ac')
    villain_hands_str = {
        "QJ suited (QhJh)": ('Qh', 'Jh'),
        "QT suited (QhTh)": ('Qh', 'Th'),
        "Q9 suited (Qh9h)": ('Qh', '9h')
    }
    
    hero_hand = [parse_card(c) for c in hero_hand_str]
    results = {}
    
    print("Calculating equity for AsAc against various red suited hands...")
    print("This will take a couple of minutes as it checks all 1,712,304 possible boards for each matchup.\n")

    for name, v_hand_str in villain_hands_str.items():
        villain_hand = [parse_card(c) for c in v_hand_str]
        equity, wins, ties, total = calculate_equity(hero_hand, villain_hand)
        results[name] = equity
        
        print(f"Matchup: {hero_hand_str[0]}{hero_hand_str[1]} vs {v_hand_str[0]}{v_hand_str[1]}")
        print(f"  - Boards simulated: {total:,}")
        print(f"  - Hero wins: {wins:,}, Hero ties: {ties:,}")
        print(f"  - Equity Equation: ({wins} + {ties} / 2) / {total}")
        print(f"  - Hero Equity: {equity:.2%}\n")

    # Find the hand that minimizes hero's equity
    min_equity_hand = min(results, key=results.get)
    
    print("--- Conclusion ---")
    print(f"Your equity is lowest against {min_equity_hand}.")
    print("Therefore, this is the hand you would least like to see.")