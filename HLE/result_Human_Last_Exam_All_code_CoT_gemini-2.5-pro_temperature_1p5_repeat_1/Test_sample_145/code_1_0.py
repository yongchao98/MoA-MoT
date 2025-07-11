import collections
import itertools
import random

def get_hand_rank(hand):
    """
    Calculates the rank of a 5-card poker hand.
    Ranks are represented as tuples for easy comparison.
    (Hand Type, Kicker1, Kicker2, ...)
    Hand Types: 8=SF, 7=Quads, 6=FH, 5=Flush, 4=Straight, 3=Trips, 2=TwoPair, 1=Pair, 0=HighCard
    """
    ranks = sorted([card[0] for card in hand], reverse=True)
    suits = [card[1] for card in hand]
    
    is_flush = len(set(suits)) == 1
    # Check for A-5 low straight
    is_ace_low_straight = (ranks == [14, 5, 4, 3, 2])
    # Check for regular straight
    is_straight = (len(set(ranks)) == 5 and ranks[0] - ranks[4] == 4) or is_ace_low_straight

    # For ranking ace-low straights, treat the Ace as low
    if is_ace_low_straight:
        eval_ranks = [5, 4, 3, 2, 1]
    else:
        eval_ranks = ranks

    if is_straight and is_flush:
        return (8, tuple(eval_ranks))  # Straight Flush
        
    counts = collections.Counter(ranks)
    sorted_counts = sorted(counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
    
    if sorted_counts[0][1] == 4: # Four of a kind
        quad_rank = sorted_counts[0][0]
        kicker = sorted_counts[1][0]
        return (7, quad_rank, kicker)

    if sorted_counts[0][1] == 3 and sorted_counts[1][1] == 2: # Full House
        trip_rank = sorted_counts[0][0]
        pair_rank = sorted_counts[1][0]
        return (6, trip_rank, pair_rank)
        
    if is_flush:
        return (5, tuple(ranks))  # Flush

    if is_straight:
        return (4, tuple(eval_ranks))  # Straight

    if sorted_counts[0][1] == 3: # Three of a kind
        trip_rank = sorted_counts[0][0]
        kickers = [r for r in ranks if r != trip_rank]
        return (3, trip_rank, kickers[0], kickers[1])
        
    if sorted_counts[0][1] == 2 and len(sorted_counts) > 1 and sorted_counts[1][1] == 2: # Two Pair
        p1 = sorted_counts[0][0]
        p2 = sorted_counts[1][0]
        kicker = [r for r in ranks if r not in [p1, p2]][0]
        return (2, max(p1, p2), min(p1, p2), kicker)

    if sorted_counts[0][1] == 2: # One Pair
        pair_rank = sorted_counts[0][0]
        kickers = [r for r in ranks if r != pair_rank]
        return (1, pair_rank, kickers[0], kickers[1], kickers[2])
        
    return (0, tuple(ranks)) # High Card

def evaluate_7_cards(seven_cards):
    """Evaluates the best 5-card hand from a set of 7 cards."""
    best_rank = (-1,)
    for five_card_combo in itertools.combinations(seven_cards, 5):
        current_rank = get_hand_rank(list(five_card_combo))
        if current_rank > best_rank:
            best_rank = current_rank
    return best_rank

def card_to_string(card):
    """Converts a card tuple to a readable string."""
    ranks = {14: 'A', 13: 'K', 12: 'Q', 11: 'J', 10: 'T'}
    rank_str = ranks.get(card[0], str(card[0]))
    suit_str = card[1]
    return f"{rank_str}{suit_str}"
    
def run_simulation(hero_hand, villain_hand, num_simulations=100000):
    """Runs a Monte Carlo simulation for a given matchup."""
    wins = 0
    ties = 0
    losses = 0
    
    # Create the deck
    ranks = range(2, 15)
    suits = ['s', 'h', 'd', 'c']
    full_deck = [(r, s) for r in ranks for s in suits]
    
    # Remove known cards from the deck
    remaining_deck = [card for card in full_deck if card not in hero_hand and card not in villain_hand]
    
    for _ in range(num_simulations):
        # Shuffle and draw a 5-card board
        random.shuffle(remaining_deck)
        board = remaining_deck[:5]
        
        # Evaluate hands
        hero_full_hand = hero_hand + board
        villain_full_hand = villain_hand + board
        
        hero_rank = evaluate_7_cards(hero_full_hand)
        villain_rank = evaluate_7_cards(villain_full_hand)
        
        # Determine outcome
        if hero_rank > villain_rank:
            wins += 1
        elif villain_rank > hero_rank:
            losses += 1
        else:
            ties += 1
            
    equity = (wins + ties / 2) / num_simulations
    
    hero_str = "".join([card_to_string(c) for c in hero_hand])
    villain_str = "".join([card_to_string(c) for c in villain_hand])
    
    print(f"Matchup: {hero_str} vs {villain_str}")
    print(f"Total simulations: {num_simulations}")
    print(f"Aces Wins: {wins}, Ties: {ties}, Losses: {losses}")
    print(f"Aces Equity = ({wins} + {ties} / 2) / {num_simulations} = {equity:.2%}\n")
    return equity

def main():
    # Define hands
    # Ranks: 14=A, 13=K, 12=Q, 11=J, 10=T, 9=9
    # Suits: s=spades, c=clubs, h=hearts, d=diamonds
    hero_hand = [(14, 's'), (14, 'c')]  # Black Aces
    
    # Red suited opponent hands
    villain_qj = [(12, 'h'), (11, 'h')] # QJ suited
    villain_qt = [(12, 'h'), (10, 'h')] # QT suited
    villain_q9 = [(12, 'h'), ( 9, 'h')] # Q9 suited
    
    villain_hands = {
      "QJ": villain_qj,
      "QT": villain_qt,
      "Q9": villain_q9
    }
    
    results = {}
    for name, hand in villain_hands.items():
        results[name] = run_simulation(hero_hand, hand)
        
    min_equity_hand = min(results, key=results.get)
    
    print("--- Conclusion ---")
    print(f"The hand that minimizes your equity as the Aces holder is {min_equity_hand} suited.")
    print("This is because it has the highest chance to win against your aces.")

if __name__ == "__main__":
    main()