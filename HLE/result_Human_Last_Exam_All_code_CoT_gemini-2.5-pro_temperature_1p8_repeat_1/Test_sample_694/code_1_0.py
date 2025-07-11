import collections
import itertools
from functools import lru_cache

@lru_cache(maxsize=None)
def get_best_score(hand_tuple):
    """
    Calculates the best possible score for a given 5-dice hand.
    Assumes all scoring categories are available.
    The hand is passed as a sorted tuple to enable memoization.
    """
    hand = list(hand_tuple)
    counts = collections.Counter(hand)
    total_sum = sum(hand)
    scores = []

    # Upper Section
    for i in range(1, 7):
        scores.append(counts[i] * i)

    # Lower Section
    # 3 and 4 of a kind
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum) # 3-of-a-kind
    else:
        scores.append(0)

    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum) # 4-of-a-kind
    else:
        scores.append(0)
    
    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    else:
        scores.append(0)

    # Straights
    unique_dice = set(hand)
    # Small Straight (4 sequential dice)
    if ({1, 2, 3, 4} <= unique_dice or
        {2, 3, 4, 5} <= unique_dice or
        {3, 4, 5, 6} <= unique_dice):
        scores.append(30)
    else:
        scores.append(0)
    
    # Large Straight (5 sequential dice)
    if unique_dice in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]:
        scores.append(40)
    else:
        scores.append(0)

    # Yahtzee
    if 5 in counts.values():
        scores.append(50)
    else:
        scores.append(0)
        
    # Chance
    scores.append(total_sum)
    
    return max(scores)

def calculate_expected_score(dice_to_keep):
    """
    Calculates the expected score for keeping a subset of dice and rerolling the rest.
    """
    num_to_reroll = 5 - len(dice_to_keep)
    
    if num_to_reroll == 0:
        return get_best_score(tuple(sorted(dice_to_keep)))

    total_score = 0
    num_outcomes = 6 ** num_to_reroll
    
    # Iterate through all possible outcomes of the dice roll
    for roll_outcome in itertools.product(range(1, 7), repeat=num_to_reroll):
        final_hand = list(dice_to_keep) + list(roll_outcome)
        total_score += get_best_score(tuple(sorted(final_hand)))
        
    return total_score / num_outcomes

def solve_yahtzee_reroll():
    """
    Main function to analyze the Yahtzee hand and determine the best move.
    """
    initial_hand = [3, 3, 3, 5, 6]

    # Generate all unique subsets of dice to keep
    # Manually defined for clarity and to handle duplicates correctly
    keep_options = [
        [],
        [3], [5], [6],
        [3, 3], [3, 5], [3, 6], [5, 6],
        [3, 3, 3], [3, 3, 5], [3, 3, 6], [3, 5, 6],
        [3, 3, 3, 5], [3, 3, 3, 6], [3, 3, 5, 6],
        [3, 3, 3, 5, 6]
    ]

    best_option = None
    max_expected_score = -1

    print("Calculating expected score for each possible 'keep' decision:")
    print("-" * 60)
    
    for option in keep_options:
        expected_score = calculate_expected_score(option)
        option_str = str(option) if option else "None (reroll all)"
        print(f"Keeping {option_str:<18}: Expected Score = {expected_score:.2f}")

        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_option = option

    print("-" * 60)
    print("\nConclusion:")
    print(f"To maximize your expected score, you should keep the dice: {best_option}")
    print(f"This yields a maximum expected score of {max_expected_score:.2f}")

if __name__ == "__main__":
    solve_yahtzee_reroll()