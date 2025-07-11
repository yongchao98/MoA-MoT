import itertools
from collections import Counter

def get_best_score(dice):
    """Calculates the best possible score for a given hand of 5 dice."""
    counts = Counter(dice)
    total_sum = sum(dice)
    scores = {}

    # Upper Section
    for i in range(1, 7):
        scores[f'upper_{i}'] = i * counts[i]

    # Lower Section
    scores['3_of_a_kind'] = total_sum if any(c >= 3 for c in counts.values()) else 0
    scores['4_of_a_kind'] = total_sum if any(c >= 4 for c in counts.values()) else 0
    scores['full_house'] = 25 if sorted(counts.values()) == [2, 3] else 0
    
    unique_dice_sorted = sorted(list(set(dice)))
    is_small_straight = False
    if len(unique_dice_sorted) >= 4:
        for i in range(len(unique_dice_sorted) - 3):
            if unique_dice_sorted[i+1] == unique_dice_sorted[i]+1 and \
               unique_dice_sorted[i+2] == unique_dice_sorted[i]+2 and \
               unique_dice_sorted[i+3] == unique_dice_sorted[i]+3:
                is_small_straight = True
                break
    scores['small_straight'] = 30 if is_small_straight else 0
    
    is_large_straight = len(unique_dice_sorted) == 5 and (unique_dice_sorted[4] - unique_dice_sorted[0] == 4)
    scores['large_straight'] = 40 if is_large_straight else 0

    scores['yahtzee'] = 50 if 5 in counts.values() else 0
    scores['chance'] = total_sum
    
    return max(scores.values())

def get_expected_score(kept_dice):
    """Calculates the average expected score for a given set of kept dice."""
    num_reroll = 5 - len(kept_dice)
    if num_reroll == 0:
        return get_best_score(kept_dice)

    total_score_sum = 0
    num_outcomes = 6 ** num_reroll
    
    possible_rolls = itertools.product(range(1, 7), repeat=num_reroll)

    for roll_outcome in possible_rolls:
        final_hand = list(kept_dice) + list(roll_outcome)
        total_score_sum += get_best_score(tuple(final_hand))

    return total_score_sum / num_outcomes

def main():
    """Analyzes the hand and prints the best strategy."""
    initial_hand = (3, 3, 3, 5, 6)

    # Define the different strategies (subsets of dice to keep).
    strategies = {
        "Keep (3, 3, 3) and reroll two dice": (3, 3, 3),
        "Keep (5, 6) and reroll three dice": (5, 6),
        "Keep (3, 5, 6) and reroll two dice": (3, 5, 6),
        "Keep all (no reroll)": (3, 3, 3, 5, 6),
        "Reroll all five dice": (),
    }
    
    results = {}
    for name, kept_dice in strategies.items():
        results[name] = get_expected_score(kept_dice)
    
    # Find and print the best strategy.
    best_strategy_name = max(results, key=results.get)
    max_score = results[best_strategy_name]
    dice_to_keep = strategies[best_strategy_name]

    print(f"Analysis for hand: {', '.join(map(str, initial_hand))}\n")
    print("Expected score for different reroll strategies:")
    print("-" * 50)
    for name, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        print(f"- {name:<38}: {score:.2f}")
    print("-" * 50)
    
    dice_to_keep_str = ', '.join(map(str, sorted(dice_to_keep))) if dice_to_keep else "nothing"
    print(f"\nTo maximize your expected score, you should keep the dice with values: {dice_to_keep_str}.")
    print(f"This strategy has the highest expected score of {max_score:.2f}.")

if __name__ == '__main__':
    main()
<<<3, 3, 3>>>