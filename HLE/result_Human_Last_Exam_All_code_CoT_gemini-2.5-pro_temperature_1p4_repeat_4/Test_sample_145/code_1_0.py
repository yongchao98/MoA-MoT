import itertools
import random

def get_card_value(card):
    """Converts a card string (e.g., 'Th', 'As') to a numerical value."""
    rank = card[:-1]
    if rank == 'T': return 10
    if rank == 'J': return 11
    if rank == 'Q': return 12
    if rank == 'K': return 13
    if rank == 'A': return 14
    return int(rank)

def evaluate_7_cards(hand_7_cards):
    """
    Evaluates the best 5-card hand from 7 cards.
    Returns a tuple representing the hand's rank for easy comparison.
    The format is (hand_type, value1, value2, ...), e.g.,
    Full House (A,A,A,K,K) -> (6, 14, 13)
    Flush (K,Q,T,5,2) -> (5, 13, 12, 10, 5, 2)
    """
    values = sorted([get_card_value(c) for c in hand_7_cards], reverse=True)
    suits = [c[-1] for c in hand_7_cards]
    
    # --- Data Pre-processing ---
    suit_counts = {}
    for card in hand_7_cards:
        suit = card[-1]
        value = get_card_value(card)
        if suit not in suit_counts:
            suit_counts[suit] = []
        suit_counts[suit].append(value)

    rank_counts = {v: values.count(v) for v in set(values)}

    # --- Check for Flush and Straight Flush ---
    flush_suit = None
    for s, s_values in suit_counts.items():
        if len(s_values) >= 5:
            flush_suit = s
            break
            
    if flush_suit:
        flush_values = sorted(suit_counts[flush_suit], reverse=True)
        unique_flush_values = sorted(list(set(flush_values)), reverse=True)
        # Check for straight within the flush cards
        # Ace-low straight check (5, 4, 3, 2, A)
        if set([14, 2, 3, 4, 5]).issubset(set(unique_flush_values)):
            return (8, (5,)) # Rank 8, high card 5
        # Regular straight check
        for i in range(len(unique_flush_values) - 4):
            if unique_flush_values[i] - unique_flush_values[i+4] == 4:
                return (8, (unique_flush_values[i],)) # Rank 8, high card of straight
        # It's a regular flush
        return (5, tuple(flush_values[:5]))
            
    # --- Check for non-flush hands (Quads, Full House, etc.) ---
    # Sort ranks by count first, then by rank value
    sorted_ranks = sorted(rank_counts.items(), key=lambda item: (-item[1], -item[0]))
    
    counts = [item[1] for item in sorted_ranks]
    ranks = [item[0] for item in sorted_ranks]
    
    # Four of a Kind
    if counts[0] == 4:
        return (7, (ranks[0], ranks[1])) # Rank 7, quad rank, kicker rank
    
    # Full House
    if counts[0] == 3 and counts[1] >= 2:
        return (6, (ranks[0], ranks[1])) # Rank 6, trips rank, pair rank

    # --- Check for Straight (must happen after flush check) ---
    unique_values = sorted(list(set(values)), reverse=True)
    # Ace-low straight check (A, 2, 3, 4, 5)
    if set([14, 2, 3, 4, 5]).issubset(set(unique_values)):
        return (4, (5,)) # Rank 4, high card 5
    # Regular straight check
    for i in range(len(unique_values) - 4):
        if unique_values[i] - unique_values[i+4] == 4:
            return (4, (unique_values[i],)) # Rank 4, high card of straight

    # Three of a Kind
    if counts[0] == 3:
        return (3, (ranks[0], ranks[1], ranks[2])) # Rank 3, trips rank, kickers

    # Two Pair
    if counts[0] == 2 and counts[1] == 2:
        return (2, (ranks[0], ranks[1], ranks[2])) # Rank 2, pairs, kicker

    # One Pair
    if counts[0] == 2:
        return (1, (ranks[0],) + tuple(ranks[1:5])) # Rank 1, pair, kickers
        
    # High Card
    return (0, tuple(ranks[:5])) # Rank 0, kickers

def calculate_equity(hero_hand, villain_hand, num_simulations=200000):
    """Calculates hero's equity against a villain via Monte Carlo simulation."""
    ranks_str = '23456789TJQKA'
    suits_str = 'shdc'
    deck = [r + s for r in ranks_str for s in suits_str]

    # Remove known cards from the deck
    for card in hero_hand + villain_hand:
        deck.remove(card)

    wins = 0
    ties = 0
    
    for _ in range(num_simulations):
        board = random.sample(deck, 5)
        
        hero_full_hand = hero_hand + board
        villain_full_hand = villain_hand + board

        hero_rank = evaluate_7_cards(hero_full_hand)
        villain_rank = evaluate_7_cards(villain_full_hand)

        if hero_rank > villain_rank:
            wins += 1
        elif hero_rank == villain_rank:
            ties += 1
            
    equity = (wins + 0.5 * ties) / num_simulations
    return equity

def main():
    """Main function to run the poker simulation and find the answer."""
    # Define the hands
    hero_hand = ['As', 'Ac']  # Two black aces
    
    # Red suited hands for the villain (we'll use hearts)
    villain_hands = {
        'QJ suited': ['Qh', 'Jh'],
        'QT suited': ['Qh', 'Th'],
        'Q9 suited': ['Qh', '9h'],
    }
    
    results = {}
    
    print("Calculating equity for Aces vs. different red suited hands...")
    print(f"Hero's hand: {hero_hand[0]} {hero_hand[1]}\n")

    for name, hand in villain_hands.items():
        equity = calculate_equity(hero_hand, hand)
        results[name] = equity
        print(f"Equity for Aces vs. {name} ({hand[0]} {hand[1]}): {equity:.2%}")

    # Find the hand that minimizes the hero's (Aces) equity
    min_equity_hand = min(results, key=results.get)
    
    print("\n--- Conclusion ---")
    print(f"The hand that minimizes the equity of the two black aces is {min_equity_hand}.")
    print("Therefore, this is the hand you would least like to see your opponent hold.")

if __name__ == "__main__":
    main()
<<<A>>>