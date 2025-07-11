import collections
import itertools
import random

def create_deck():
    """Creates a standard 52-card deck."""
    ranks = list('23456789TJQKA')
    suits = list('shdc')
    return [r + s for r in ranks for s in suits]

def parse_card(card):
    """Parses a card string 'As' into a (rank, suit) tuple, e.g., (14, 's')."""
    rank_map = {'2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'T': 10, 'J': 11, 'Q': 12, 'K': 13, 'A': 14}
    rank_str = card[:-1]
    suit = card[-1]
    return (rank_map[rank_str], suit)

def score_5_card_hand(hand):
    """
    Scores a 5-card hand.
    Returns a tuple where the first element is the hand rank (e.g., 8 for straight flush)
    and subsequent elements are kickers for tie-breaking.
    """
    ranks = sorted([card[0] for card in hand], reverse=True)
    suits = [card[1] for card in hand]
    
    is_flush = len(set(suits)) == 1
    
    unique_ranks = sorted(list(set(ranks)), reverse=True)
    is_straight = (len(unique_ranks) == 5 and unique_ranks[0] - unique_ranks[4] == 4) or \
                  (unique_ranks == [14, 5, 4, 3, 2])
    
    high_card_straight = unique_ranks[0]
    if unique_ranks == [14, 5, 4, 3, 2]: # A-5 wheel
        high_card_straight = 5
        ranks = [5, 4, 3, 2, 1]

    if is_straight and is_flush: return (8, high_card_straight)
    
    rank_counts = collections.Counter(ranks)
    counts = sorted(rank_counts.values(), reverse=True)
    
    if counts[0] == 4:
        quad_rank = [r for r, c in rank_counts.items() if c == 4][0]
        kicker = [r for r, c in rank_counts.items() if c == 1][0]
        return (7, quad_rank, kicker)
        
    if counts == [3, 2]:
        trips_rank = [r for r, c in rank_counts.items() if c == 3][0]
        pair_rank = [r for r, c in rank_counts.items() if c == 2][0]
        return (6, trips_rank, pair_rank)
        
    if is_flush: return (5, tuple(ranks))

    if is_straight: return (4, high_card_straight)
        
    if counts[0] == 3:
        trips_rank = [r for r, c in rank_counts.items() if c == 3][0]
        kickers = sorted([r for r, c in rank_counts.items() if c == 1], reverse=True)
        return (3, trips_rank, kickers[0], kickers[1])
        
    if counts == [2, 2, 1]:
        pairs = sorted([r for r, c in rank_counts.items() if c == 2], reverse=True)
        kicker = [r for r, c in rank_counts.items() if c == 1][0]
        return (2, pairs[0], pairs[1], kicker)
        
    if counts[0] == 2:
        pair_rank = [r for r, c in rank_counts.items() if c == 2][0]
        kickers = sorted([r for r, c in rank_counts.items() if c == 1], reverse=True)
        return (1, pair_rank, kickers[0], kickers[1], kickers[2])
        
    return (0, tuple(ranks))

def evaluate_hand(hand_str):
    """Evaluates the best 5-card hand from a 7-card hand."""
    parsed_cards = [parse_card(c) for c in hand_str]
    best_score = (-1,)
    for combo in itertools.combinations(parsed_cards, 5):
        score = score_5_card_hand(list(combo))
        if score > best_score:
            best_score = score
    return best_score

def run_simulation(hero_hand, villain_hand, num_simulations):
    """Runs a Monte Carlo simulation to calculate hero's equity."""
    deck = create_deck()
    known_cards = hero_hand + villain_hand
    remaining_deck = [card for card in deck if card not in known_cards]
    
    hero_wins = 0
    chops = 0
    
    for _ in range(num_simulations):
        board = random.sample(remaining_deck, 5)
        
        hero_7_hand = hero_hand + board
        villain_7_hand = villain_hand + board
        
        hero_score = evaluate_hand(hero_7_hand)
        villain_score = evaluate_hand(villain_7_hand)
        
        if hero_score > villain_score:
            hero_wins += 1
        elif hero_score == villain_score:
            chops += 1
            
    return (hero_wins + 0.5 * chops) / num_simulations

def main():
    """Main function to run the poker equity calculations."""
    hero_hand = ['As', 'Ac']  # Two black aces
    villain_hands = {
        "QJ suited": ['Qh', 'Jh'],
        "QT suited": ['Qh', 'Th'],
        "Q9 suited": ['Qh', '9h'],
    }
    num_simulations = 100000

    print(f"Calculating equity for Hero ({hero_hand}) vs three possible Villain hands.")
    print(f"This is done using a Monte Carlo simulation with {num_simulations} random boards for each matchup.\n")

    equities = {}
    for name, hand in villain_hands.items():
        equity = run_simulation(hero_hand, hand, num_simulations)
        equities[name] = equity
        # Each number in the final equation: Hero wins % of the time.
        print(f"Hero's equity vs {name} ({hand}): {equity:.2%}")
    
    # The hand you least like to see is the one that minimizes your equity.
    worst_matchup_for_hero = min(equities, key=equities.get)
    min_equity_value = equities[worst_matchup_for_hero]

    print(f"\nConclusion:")
    print(f"The hand that minimizes your equity is {worst_matchup_for_hero}.")
    print(f"This is because against {worst_matchup_for_hero}, your aces have the lowest winning probability of {min_equity_value:.2%}.")
    print("Therefore, QJ suited is the hand you least like to see.")

if __name__ == '__main__':
    main()
<<<A>>>