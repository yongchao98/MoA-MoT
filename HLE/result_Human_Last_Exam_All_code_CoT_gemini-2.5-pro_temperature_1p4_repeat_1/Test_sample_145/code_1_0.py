import itertools
from collections import Counter

def get_hand_rank(seven_cards):
    """
    Evaluates 7 cards and returns the best 5-card hand rank as a tuple.
    The tuple format is (rank_value, kicker_values...) which allows for direct comparison.
    Higher rank_value is better.
    Rank values: 8=StraightFlush, 7=Quads, 6=FullHouse, 5=Flush, 4=Straight, 
                   3=Trips, 2=TwoPair, 1=Pair, 0=HighCard.
    """
    # Card ranks mapping '2'->2, ..., 'A'->14
    rank_map = {r: i for i, r in enumerate('23456789TJQKA', 2)}
    
    ranks = [rank_map[c[0]] for c in seven_cards]
    suits = [c[1] for c in seven_cards]

    # --- Hand Checks ---
    # Check for Flush and Straight Flush
    suit_counts = Counter(suits)
    flush_suit = None
    if suit_counts:
        for suit, count in suit_counts.items():
            if count >= 5:
                flush_suit = suit
                break

    is_flush = flush_suit is not None
    is_straight = False
    straight_top_rank = 0

    # Get unique sorted ranks for straight check
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    # Ace can be low for A-2-3-4-5 straight
    if 14 in unique_ranks:
        unique_ranks.append(1)

    # Check for straight
    for i in range(len(unique_ranks) - 4):
        # Check for 5 consecutive ranks
        if unique_ranks[i] - 4 == unique_ranks[i+4]:
            is_straight = True
            # The top card of A-2-3-4-5 is 5, not the Ace
            straight_top_rank = unique_ranks[i] if unique_ranks[i] > 5 else 5
            break
    
    # Check Straight Flush
    if is_flush and is_straight:
        flush_ranks = sorted([rank_map[c[0]] for c in seven_cards if c[1] == flush_suit], reverse=True)
        # Check for straight within the flush cards
        if 14 in flush_ranks:
            flush_ranks.append(1)
            
        for i in range(len(flush_ranks) - 4):
            if flush_ranks[i] - 4 == flush_ranks[i+4]:
                sf_top_rank = flush_ranks[i] if flush_ranks[i] > 5 else 5
                return (8, sf_top_rank) # Straight Flush

    # Count rank occurrences for Quads, Full House, Trips, Pairs
    rank_counts = Counter(ranks)
    # Sort by count (desc), then rank (desc)
    sorted_counts = sorted(rank_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)

    # Check for Quads
    if sorted_counts[0][1] == 4:
        quad_rank = sorted_counts[0][0]
        kickers = sorted([r for r in ranks if r != quad_rank], reverse=True)
        return (7, quad_rank, kickers[0]) # Four of a Kind
        
    # Check for Full House
    if sorted_counts[0][1] == 3 and sorted_counts[1][1] >= 2:
        return (6, sorted_counts[0][0], sorted_counts[1][0]) # Full House

    # Check for Flush
    if is_flush:
        flush_ranks = sorted([rank_map[c[0]] for c in seven_cards if c[1] == flush_suit], reverse=True)
        return (5, *flush_ranks[:5]) # Flush

    # Check for Straight
    if is_straight:
        return (4, straight_top_rank) # Straight

    # Check for Three of a Kind
    if sorted_counts[0][1] == 3:
        trip_rank = sorted_counts[0][0]
        kickers = sorted([r for r in ranks if r != trip_rank], reverse=True)
        return (3, trip_rank, *kickers[:2]) # Three of a Kind

    # Check for Two Pair
    if sorted_counts[0][1] == 2 and sorted_counts[1][1] == 2:
        pair_ranks = sorted([sorted_counts[0][0], sorted_counts[1][0]], reverse=True)
        kicker_rank = sorted([r for r in ranks if r not in pair_ranks], reverse=True)[0]
        return (2, *pair_ranks, kicker_rank) # Two Pair

    # Check for One Pair
    if sorted_counts[0][1] == 2:
        pair_rank = sorted_counts[0][0]
        kickers = sorted([r for r in ranks if r != pair_rank], reverse=True)
        return (1, pair_rank, *kickers[:3]) # One Pair

    # High Card
    return (0, *sorted(ranks, reverse=True)[:5])


def calculate_equity(hero_hand, villain_hand):
    """Calculates the equity for the hero's hand against the villain's."""
    # Create a standard 52-card deck
    ranks = '23456789TJQKA'
    suits = 'shdc' # spades, hearts, diamonds, clubs
    deck = {r + s for r in ranks for s in suits}

    # Remove known cards from the deck
    deck.remove(hero_hand[0])
    deck.remove(hero_hand[1])
    deck.remove(villain_hand[0])
    deck.remove(villain_hand[1])
    
    wins = 0
    ties = 0
    total_hands = 0

    # Iterate over all possible 5-card boards
    for board in itertools.combinations(deck, 5):
        total_hands += 1
        hero_seven_cards = hero_hand + board
        villain_seven_cards = villain_hand + board
        
        hero_rank = get_hand_rank(hero_seven_cards)
        villain_rank = get_hand_rank(villain_seven_cards)
        
        if hero_rank > villain_rank:
            wins += 1
        elif hero_rank == villain_rank:
            ties += 1

    equity = (wins + 0.5 * ties) / total_hands
    return wins, ties, total_hands, equity


def main():
    """Main function to run the simulation and print results."""
    # Hero hand: Two black aces (Spades and Clubs)
    hero_hand = ('As', 'Ac')
    
    # Villain hands: Suited red cards (Diamonds)
    villain_hands_to_test = {
        "QJ suited": ('Qd', 'Jd'),
        "QT suited": ('Qd', 'Td'),
        "Q9 suited": ('Qd', '9d'),
    }

    results = {}
    
    print("Calculating equity for Black Aces vs. various suited red hands...\n")

    for name, villain_hand in villain_hands_to_test.items():
        print(f"--- Calculating for {hero_hand[0]}{hero_hand[1]} vs {villain_hand[0]}{villain_hand[1]} ---")
        wins, ties, total, equity = calculate_equity(hero_hand, villain_hand)
        
        results[name] = equity
        
        print(f"Total boards simulated: {total}")
        print(f"Hero wins: {wins}")
        print(f"Ties: {ties}")
        print(f"Hero losses: {total - wins - ties}")
        print(f"Equity Equation: (Wins + 0.5 * Ties) / Total")
        print(f"Equity = ({wins} + 0.5 * {ties}) / {total} = {equity:.4f} or {equity:.2%}")
        print("-" * 30 + "\n")

    # Find the hand that minimizes the hero's (Aces) equity
    min_equity_hand = min(results, key=results.get)
    min_equity_value = results[min_equity_hand]

    print("\n--- Conclusion ---")
    print("The equity for Aces is lowest against the hand that is most connected,")
    print("as it has the best chance of making a straight to beat the aces.")
    print(f"The hand that minimizes your equity is: {min_equity_hand}")
    print(f"Your equity against {min_equity_hand} is only {min_equity_value:.2%}.")

if __name__ == '__main__':
    main()
<<<A>>>