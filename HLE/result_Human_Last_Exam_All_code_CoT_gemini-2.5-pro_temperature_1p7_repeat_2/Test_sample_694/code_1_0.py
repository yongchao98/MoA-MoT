import itertools
from collections import Counter

def calculate_best_score(hand):
    """
    Calculates the maximum possible score for a given 5-dice hand
    by checking all scoring categories.
    """
    if not isinstance(hand, list) or len(hand) != 5:
        return 0
        
    counts = Counter(hand)
    total_sum = sum(hand)
    scores = {}

    # Upper Section Scores
    for i in range(1, 7):
        scores[f'upper_{i}'] = counts.get(i, 0) * i

    # Chance
    scores['chance'] = total_sum

    # Yahtzee (50 points)
    if 5 in counts.values():
        scores['yahtzee'] = 50

    # Large Straight (40 points)
    unique_dice = sorted(counts.keys())
    if unique_dice in [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6]]:
        scores['large_straight'] = 40

    # Small Straight (30 points)
    is_small_straight = False
    if len(unique_dice) >= 4:
        for i in range(len(unique_dice) - 3):
            if unique_dice[i+1] == unique_dice[i] + 1 and \
               unique_dice[i+2] == unique_dice[i] + 2 and \
               unique_dice[i+3] == unique_dice[i] + 3:
                is_small_straight = True
                break
    if is_small_straight:
        scores['small_straight'] = 30

    # Full House (25 points)
    if sorted(counts.values()) == [2, 3]:
        scores['full_house'] = 25

    # Four of a Kind (sum of all dice)
    if any(c >= 4 for c in counts.values()):
        scores['four_of_a_kind'] = total_sum

    # Three of a Kind (sum of all dice)
    if any(c >= 3 for c in counts.values()):
        scores['three_of_a_kind'] = total_sum

    return max(scores.values()) if scores else 0

def calculate_expected_score(dice_to_keep):
    """
    Calculates the expected score for keeping a given set of dice and rerolling the rest.
    """
    num_to_reroll = 5 - len(dice_to_keep)
    
    if num_to_reroll == 0:
        return calculate_best_score(dice_to_keep)

    total_score = 0
    num_outcomes = 6 ** num_to_reroll
    
    # Generate all possible outcomes for the dice being rerolled
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_to_reroll):
        final_hand = dice_to_keep + list(reroll_outcome)
        total_score += calculate_best_score(final_hand)
        
    return total_score / num_outcomes

def main():
    """
    Main function to analyze the Yahtzee hand and suggest the best strategy.
    """
    initial_hand = [3, 3, 3, 5, 6]
    
    strategies = {
        "Keep all dice [3, 3, 3, 5, 6]": initial_hand,
        "Keep three-of-a-kind [3, 3, 3]": [3, 3, 3],
        "Keep highest die [6]": [6]
    }
    
    results = {}
    print(f"Analyzing hand: {initial_hand}\n")

    for desc, keep_hand in strategies.items():
        expected_score = calculate_expected_score(keep_hand)
        results[desc] = expected_score
        print(f"Strategy: {desc}")
        print(f" -> Expected Score: {expected_score:.2f}\n")
    
    best_strategy = max(results, key=results.get)
    best_score = results[best_strategy]
    
    # Extract the numbers from the best strategy description
    best_dice_to_keep = "".join(filter(lambda x: x if x in '0123456789,' else ' ', best_strategy.split('[')[1]))
    
    print("--- Conclusion ---")
    print(f"Comparing the expected scores: {list(results.values())}")
    print(f"The optimal strategy is '{best_strategy}' with an expected score of {best_score:.2f}.")
    print(f"Therefore, you should keep the dice: {best_dice_to_keep}")


if __name__ == "__main__":
    main()