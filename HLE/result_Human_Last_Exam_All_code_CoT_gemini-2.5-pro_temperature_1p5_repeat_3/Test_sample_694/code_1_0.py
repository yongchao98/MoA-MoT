import collections
import itertools
from math import pow

def get_best_score(hand):
    """
    Calculates the maximum possible score for a given 5-dice hand,
    assuming all scoring categories are open.
    """
    counts = collections.Counter(hand)
    hand_sum = sum(hand)
    best_score = 0

    # Upper Section Score
    for i in range(1, 7):
        best_score = max(best_score, i * counts.get(i, 0))

    # Lower Section Scores
    is_three_of_kind = any(c >= 3 for c in counts.values())
    is_four_of_kind = any(c >= 4 for c in counts.values())
    is_yahtzee = any(c == 5 for c in counts.values())
    is_full_house = sorted(counts.values()) == [2, 3]

    if is_three_of_kind:
        best_score = max(best_score, hand_sum)
    if is_four_of_kind:
        best_score = max(best_score, hand_sum)
    if is_full_house:
        best_score = max(best_score, 25)
    if is_yahtzee:
        best_score = max(best_score, 50)

    # Straights
    unique_dice_str = "".join(map(str, sorted(counts.keys())))
    if "1234" in unique_dice_str or "2345" in unique_dice_str or "3456" in unique_dice_str:
        best_score = max(best_score, 30)
    if len(counts) == 5 and (max(counts.keys()) - min(counts.keys()) == 4):
         best_score = max(best_score, 40)

    # Chance
    best_score = max(best_score, hand_sum)

    return best_score

memo = {}

def get_expected_score(kept_dice):
    """
    Calculates the expected score for a given set of kept dice by
    simulating all possible outcomes of the rerolled dice.
    """
    kept_dice_tuple = tuple(sorted(kept_dice))
    num_to_roll = 5 - len(kept_dice_tuple)

    if num_to_roll == 0:
        return get_best_score(kept_dice_tuple)

    if kept_dice_tuple in memo:
        return memo[kept_dice_tuple]

    total_score = 0
    num_outcomes = pow(6, num_to_roll)
    
    # Iterate through all possible outcomes for the dice being rolled
    for rolls in itertools.product(range(1, 7), repeat=num_to_roll):
        final_hand = list(kept_dice_tuple) + list(rolls)
        total_score += get_best_score(final_hand)

    expected_value = total_score / num_outcomes
    memo[kept_dice_tuple] = expected_value
    return expected_value

def main():
    """
    Analyzes the Yahtzee hand and determines the optimal play.
    """
    initial_hand = [3, 3, 3, 5, 6]

    # Generate all unique subsets of dice to keep
    subsets_to_check = set()
    for i in range(len(initial_hand) + 1):
        for subset in itertools.combinations(initial_hand, i):
            subsets_to_check.add(tuple(sorted(subset)))

    # Calculate the expected score for each possible play
    results = {}
    for subset in subsets_to_check:
        expected_score = get_expected_score(list(subset))
        results[subset] = expected_score

    # Find the best strategy
    best_strategy_dice = max(results, key=results.get)
    max_score = results[best_strategy_dice]

    # Sort and print the top strategies for context
    print("--- Yahtzee Strategy Analysis ---")
    print(f"Initial Hand: {initial_hand}\n")
    print("Top 10 Optimal Plays (sorted by expected score):")
    print("-" * 55)
    print(f"{'Dice to Keep':<25} | {'Dice to Reroll':<15} | {'Expected Score':>15}")
    print("-" * 55)

    sorted_results = sorted(results.items(), key=lambda item: item[1], reverse=True)
    for kept, ev in sorted_results[:10]:
        kept_str = "Nothing" if not kept else str(list(kept))
        reroll_count = 5 - len(kept)
        print(f"{kept_str:<25} | {reroll_count:<15} | {ev:>15.2f}")
    
    print("-" * 55)

    # Print the final recommendation
    final_answer_dice = list(best_strategy_dice)
    if not final_answer_dice:
      final_answer_str = "nothing and reroll all five dice"
    else:
      final_answer_str = f"the dice with values: {final_answer_dice}"

    print(f"\nConclusion:")
    print(f"To maximize your expected score, you should keep {final_answer_str}.")
    print(f"This strategy has the highest expected score of {max_score:.2f} points.")
    
    # Format the final answer string for the tag
    final_answer_tag_content = ', '.join(map(str, sorted(final_answer_dice)))
    print(f"\n<<<[{final_answer_tag_content}]>>>")


if __name__ == '__main__':
    main()