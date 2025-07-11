import collections
import itertools

def get_best_score(dice):
    """Calculates the highest possible score for a given 5-dice hand."""
    counts = collections.Counter(dice)
    total_sum = sum(dice)
    scores = []

    # Upper Section
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Three-Of-A-Kind: Sum of all dice
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Four-Of-A-Kind: Sum of all dice
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Full House: 25 points
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    else:
        scores.append(0)

    # Small Straight: 30 points (sequence of 4)
    unique_dice_set = set(dice)
    small_straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    is_small_straight = any(s.issubset(unique_dice_set) for s in small_straights)
    scores.append(30 if is_small_straight else 0)

    # Large Straight: 40 points (sequence of 5)
    is_large_straight = len(unique_dice_set) == 5 and (max(unique_dice_set) - min(unique_dice_set) == 4)
    scores.append(40 if is_large_straight else 0)

    # Yahtzee: 50 points
    if 5 in counts.values():
        scores.append(50)
    else:
        scores.append(0)

    # Chance: Sum of all dice
    scores.append(total_sum)

    return max(scores)

def calculate_expected_score(kept_dice):
    """Calculates the expected score for a given set of kept dice."""
    num_reroll = 5 - len(kept_dice)
    if num_reroll == 0:
        return get_best_score(kept_dice)

    total_score = 0
    num_outcomes = 6 ** num_reroll
    
    # Iterate through all possible outcomes of the reroll
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
        final_hand = kept_dice + list(reroll_outcome)
        total_score += get_best_score(final_hand)

    return total_score / num_outcomes

def solve_yahtzee_reroll():
    """Main function to analyze the hand and determine the best move."""
    initial_hand = [3, 3, 3, 5, 6]

    # Generate all unique subsets of dice to keep
    keep_options = set()
    for r in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, r):
            keep_options.add(tuple(sorted(combo)))

    # Calculate expected score for each option
    results = {}
    for keep_tuple in keep_options:
        expected_score = calculate_expected_score(list(keep_tuple))
        results[keep_tuple] = expected_score
    
    # Sort results for display, from best to worst
    sorted_results = sorted(results.items(), key=lambda item: item[1], reverse=True)

    print(f"Analysis for hand: {initial_hand}\n")
    print("Expected score for each possible keep strategy:")
    for keep, score in sorted_results:
        # Each number is printed in the final evaluation of strategies
        print(f"  Keep {str(list(keep)):<15}: Expected Score = {score:.2f}")

    best_keep_tuple = sorted_results[0][0]
    best_score = sorted_results[0][1]
    best_keep_list = list(best_keep_tuple)

    print("\n----------------------------------------------------")
    print("Conclusion:")
    print(f"To maximize your expected score, you should keep the dice: {best_keep_list}")
    print(f"This strategy has an expected score of {best_score:.2f} points.")
    print("----------------------------------------------------\n")

    # Final answer format
    print(f'<<<{best_keep_list}>>>')

# Run the solver
solve_yahtzee_reroll()