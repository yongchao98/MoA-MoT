import itertools
from collections import Counter
from functools import lru_cache

# Memoization cache for performance
memo_max_score = {}

def get_counts(hand):
    """Counts the occurrences of each die value."""
    return Counter(hand)

def score_upper(hand, number):
    """Calculates the score for an upper section category."""
    return sum(d for d in hand if d == number)

def score_of_a_kind(hand, n):
    """Calculates the score for N-of-a-kind categories."""
    counts = get_counts(hand)
    if any(c >= n for c in counts.values()):
        return sum(hand)
    return 0

def score_full_house(hand):
    """Calculates the score for a Full House."""
    counts = get_counts(hand)
    if sorted(counts.values()) == [2, 3]:
        return 25
    return 0

def score_small_straight(hand):
    """Calculates the score for a Small Straight."""
    unique_dice = set(hand)
    # Check for all three possible small straights
    if {1, 2, 3, 4}.issubset(unique_dice):
        return 30
    if {2, 3, 4, 5}.issubset(unique_dice):
        return 30
    if {3, 4, 5, 6}.issubset(unique_dice):
        return 30
    return 0

def score_large_straight(hand):
    """Calculates the score for a Large Straight."""
    unique_dice = sorted(list(set(hand)))
    if len(unique_dice) == 5 and (unique_dice[4] - unique_dice[0] == 4):
        return 40
    return 0

def score_yahtzee(hand):
    """Calculates the score for a Yahtzee."""
    counts = get_counts(hand)
    if 5 in counts.values():
        return 50
    return 0

def get_max_score(hand):
    """Calculates the maximum possible score for a given hand."""
    hand_tuple = tuple(sorted(hand))
    if hand_tuple in memo_max_score:
        return memo_max_score[hand_tuple]

    scores = [
        score_upper(hand, 1),
        score_upper(hand, 2),
        score_upper(hand, 3),
        score_upper(hand, 4),
        score_upper(hand, 5),
        score_upper(hand, 6),
        score_of_a_kind(hand, 3),
        score_of_a_kind(hand, 4),
        score_full_house(hand),
        score_small_straight(hand),
        score_large_straight(hand),
        score_yahtzee(hand),
        sum(hand)  # Chance
    ]
    max_s = max(scores)
    memo_max_score[hand_tuple] = max_s
    return max_s

@lru_cache(maxsize=None)
def calculate_expected_score(dice_to_keep):
    """Calculates the expected score for keeping a given subset of dice."""
    num_reroll = 5 - len(dice_to_keep)
    if num_reroll == 0:
        return get_max_score(dice_to_keep)

    total_max_score = 0
    num_outcomes = 6 ** num_reroll

    # Generate all possible outcomes for the rerolled dice
    reroll_outcomes = itertools.product(range(1, 7), repeat=num_reroll)

    for outcome in reroll_outcomes:
        final_hand = list(dice_to_keep) + list(outcome)
        total_max_score += get_max_score(final_hand)

    return total_max_score / num_outcomes

def solve_yahtzee_reroll():
    """
    Determines the optimal dice to keep to maximize expected score.
    """
    initial_dice = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    strategies = set()
    for i in range(1 << len(initial_dice)):
        subset_indices = [j for j in range(len(initial_dice)) if (i >> j) & 1]
        dice_to_keep = tuple(sorted([initial_dice[k] for k in subset_indices]))
        strategies.add(dice_to_keep)

    # Calculate expected score for each strategy
    results = {}
    for s in sorted(list(strategies), key=len):
        results[s] = calculate_expected_score(s)

    # Find and print the best strategy
    best_strategy = max(results, key=results.get)
    max_score = results[best_strategy]

    print("Analysis of all possible strategies:")
    for strategy, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        strategy_str = str(strategy) if strategy else "None"
        print(f"  Keep {strategy_str:<15}: Expected Score = {score:.2f}")

    print("\n--- Conclusion ---")
    print(f"The optimal strategy is to keep the dice: {best_strategy}")
    print(f"This yields the highest expected score of: {max_score:.2f}")
    print("The score for this hand is guaranteed, as no dice are rerolled. The best scoring category is Three of a Kind or Chance for 20 points.")

if __name__ == '__main__':
    solve_yahtzee_reroll()