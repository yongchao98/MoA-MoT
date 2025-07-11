import itertools
from collections import Counter

def score_three_of_a_kind(hand):
    """Calculates the score for Three of a Kind."""
    counts = Counter(hand)
    if any(c >= 3 for c in counts.values()):
        return sum(hand)
    return 0

def score_four_of_a_kind(hand):
    """Calculates the score for Four of a Kind."""
    counts = Counter(hand)
    if any(c >= 4 for c in counts.values()):
        return sum(hand)
    return 0

def score_full_house(hand):
    """Calculates the score for a Full House."""
    counts = Counter(hand)
    return 25 if sorted(counts.values()) == [2, 3] else 0

def score_small_straight(hand):
    """Calculates the score for a Small Straight (sequence of 4)."""
    unique_dice = set(hand)
    if {1,2,3,4}.issubset(unique_dice) or \
       {2,3,4,5}.issubset(unique_dice) or \
       {3,4,5,6}.issubset(unique_dice):
        return 30
    return 0

def score_large_straight(hand):
    """Calculates the score for a Large Straight (sequence of 5)."""
    return 40 if len(set(hand)) == 5 and max(hand) - min(hand) == 4 else 0

def score_yahtzee(hand):
    """Calculates the score for a Yahtzee (five of a kind)."""
    return 50 if len(set(hand)) == 1 else 0

def calculate_best_score(hand):
    """Calculates the maximum possible score for a given hand."""
    scores = [
        # Upper section scores
        hand.count(1) * 1,
        hand.count(2) * 2,
        hand.count(3) * 3,
        hand.count(4) * 4,
        hand.count(5) * 5,
        hand.count(6) * 6,
        # Lower section scores
        score_three_of_a_kind(hand),
        score_four_of_a_kind(hand),
        score_full_house(hand),
        score_small_straight(hand),
        score_large_straight(hand),
        score_yahtzee(hand),
        # Chance
        sum(hand)
    ]
    return max(scores)

def calculate_expected_score(kept_dice_tuple):
    """Calculates the expected score for keeping a given subset of dice."""
    kept_dice = list(kept_dice_tuple)
    num_reroll = 5 - len(kept_dice)

    if num_reroll == 0:
        return calculate_best_score(kept_dice)

    total_score = 0
    num_outcomes = 6 ** num_reroll

    # Generate all possible outcomes for the dice being rerolled
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
        final_hand = kept_dice + list(reroll_outcome)
        total_score += calculate_best_score(final_hand)

    return total_score / num_outcomes

def solve_yahtzee_reroll():
    """Main function to analyze the hand and determine the best strategy."""
    initial_hand = [3, 3, 3, 5, 6]
    possible_keeps = set()

    # Generate all unique subsets of the hand to consider keeping
    for i in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, i):
            possible_keeps.add(tuple(sorted(combo)))

    # Calculate expected score for each strategy
    results = {}
    for keep_candidate in possible_keeps:
        expected_score = calculate_expected_score(keep_candidate)
        results[keep_candidate] = expected_score

    # Sort the results to find the best strategies
    sorted_results = sorted(results.items(), key=lambda item: item[1], reverse=True)

    print(f"Analysis for the hand: {initial_hand}\n")
    print("The top 5 strategies based on expected score are:\n")

    for i in range(min(5, len(sorted_results))):
        keep, score = sorted_results[i]
        # Format the 'keep' tuple for printing
        keep_str = str(list(keep)) if keep else "None (reroll all)"
        print(f"Keep: {keep_str:<17} | Expected Score: {score:.2f}")

    best_strategy_dice = sorted_results[0][0]
    best_strategy_score = sorted_results[0][1]
    
    # Print the final equation for the best strategy
    final_dice_to_keep = ', '.join(map(str, best_strategy_dice))
    print("\n---")
    print("Final Conclusion:")
    print(f"To maximize your score, you should keep the dice: {final_dice_to_keep}")
    print(f"This strategy has the highest expected score of {best_strategy_score:.2f} points.")
    # The 'final equation' is the list of numbers to keep for the optimal outcome.
    print(f"Final numbers to keep: {final_dice_to_keep}")

# Run the analysis
solve_yahtzee_reroll()