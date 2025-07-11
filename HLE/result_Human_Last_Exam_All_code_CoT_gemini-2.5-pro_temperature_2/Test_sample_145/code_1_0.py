import itertools
from collections import Counter

def create_deck():
    """Creates a standard 52-card deck."""
    ranks = '23456789TJQKA'
    suits = 'shdc'  # spades, hearts, diamonds, clubs
    return [r + s for r in ranks for s in suits]

def rank_to_value(card_rank):
    """Converts a card rank character to its numerical value."""
    if card_rank.isdigit():
        return int(card_rank)
    return {'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}[card_rank]

def evaluate_five_card_hand(hand):
    """
    Evaluates a 5-card hand and returns a tuple representing its rank.
    Higher tuples represent better hands.
    Example Ranks: (9, [14,13,12,11,10]) for Royal Flush
                   (8, 10, 14) for Four of a Kind 10s with Ace kicker
                   (1, [14,12,10,5,3]) for Ace-high
    """
    ranks = sorted([rank_to_value(c[0]) for c in hand], reverse=True)
    suits = [c[1] for c in hand]
    is_flush = len(set(suits)) == 1

    # Check for A-5 "wheel" straight special case
    is_straight = False
    if ranks == [14, 5, 4, 3, 2]:
        is_straight = True
        rank_values = [5, 4, 3, 2, 1]  # Treat Ace as 1 for rank comparison
    else:
        unique_ranks = sorted(list(set(ranks)), reverse=True)
        if len(unique_ranks) >= 5:
           for i in range(len(unique_ranks) - 4):
               if unique_ranks[i] - unique_ranks[i+4] == 4:
                   is_straight = True
                   rank_values = unique_ranks[i:i+5]
                   break
    
    # Generate rank tuple for comparison
    if is_straight and is_flush: return (9, rank_values) # Straight Flush
    
    # Get counts of each rank for quads, full house, etc.
    rank_counts = Counter(ranks)
    sorted_counts = sorted(rank_counts.items(), key=lambda item: (item[1], item[0]), reverse=True)
    
    vals = [item[0] for item in sorted_counts]
    counts = [item[1] for item in sorted_counts]

    if counts[0] == 4: return (8, vals[0], vals[1])  # Quads
    if counts == [3, 2]: return (7, vals[0], vals[1]) # Full House
    if is_flush: return (6, ranks) # Flush
    if is_straight: return (5, rank_values) # Straight
    if counts[0] == 3: return (4, vals[0], vals[1], vals[2]) # Three of a Kind
    if counts == [2, 2, 1]: return (3, vals[0], vals[1], vals[2]) # Two Pair
    if counts[0] == 2: return (2, vals[0], vals[1], vals[2], vals[3]) # One Pair
    return (1, ranks) # High Card

def find_best_hand(hole_cards, board):
    """Takes 2 hole cards and 5 board cards and finds the best 5-card hand."""
    all_seven_cards = hole_cards + list(board)
    best_rank = (0,)
    for five_card_combo in itertools.combinations(all_seven_cards, 5):
        current_rank = evaluate_five_card_hand(list(five_card_combo))
        if current_rank > best_rank:
            best_rank = current_rank
    return best_rank

def run_simulation():
    """Main function to run the simulation and print results."""
    hero_hand = ['As', 'Ac']  # Two black aces
    villain_hands = {
        'QJ': ['Qh', 'Jh'],
        'QT': ['Qh', 'Th'],
        'Q9': ['Qh', '9h'],
    }
    
    print("Calculating exact equity of AsAc vs three different suited hands.")
    print("This will iterate through all 1,712,304 possible boards and may take a couple of minutes.")
    
    results = {}
    for name, villain_hand in villain_hands.items():
        deck = create_deck()
        hand_cards = hero_hand + villain_hand
        for card in hand_cards:
            deck.remove(card)
        
        wins, ties, total_boards = 0, 0, 0
        
        print(f"\nSimulating Hero ({''.join(hero_hand)}) vs. Villain ({''.join(villain_hand)})...")
        for board in itertools.combinations(deck, 5):
            total_boards += 1
            if total_boards % 500000 == 0:
                print(f"  ...processed {total_boards} boards...")
            
            hero_rank = find_best_hand(hero_hand, board)
            villain_rank = find_best_hand(villain_hand, board)
            
            if hero_rank > villain_rank:
                wins += 1
            elif hero_rank == villain_rank:
                ties += 1
        
        equity = (wins + (ties / 2)) / total_boards
        results[name] = equity
        print("Calculation complete.")
        print("Final Equation:")
        print(f"({wins} wins + {ties} ties / 2) / {total_boards} total boards = {equity:.4%} equity")

    min_equity_val = min(results.values())
    least_liked_hand = [name for name, equity in results.items() if equity == min_equity_val][0]

    print("\n--- Final Answer ---")
    print(f"The equities for Hero (AsAc) are:")
    for hand, equity in results.items():
        print(f"  vs. {hand}: {equity:.4%}")
    
    print(f"\nThe hand you least like to see is {least_liked_hand} because it minimizes your equity.")

if __name__ == '__main__':
    run_simulation()