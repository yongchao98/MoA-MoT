import itertools
from collections import Counter

def calculate_best_score(hand):
    """
    Calculates the best possible score for a given 5-dice hand in Yahtzee.
    This function considers all scoring categories and returns the maximum possible score for the hand.
    Note: The 35-point upper section bonus is a whole-game feature and is not included
    in single-turn expected value calculations.
    """
    scores = []
    counts = Counter(hand)
    total_sum = sum(hand)
    
    # Upper Section Scores
    for i in range(1, 7):
        scores.append(i * counts.get(i, 0))

    # Three of a Kind
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Four of a Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    else:
        scores.append(0)

    # Small Straight (sequence of 4)
    unique_dice = set(hand)
    if {1,2,3,4} <= unique_dice or {2,3,4,5} <= unique_dice or {3,4,5,6} <= unique_dice:
        scores.append(30)
    else:
        scores.append(0)

    # Large Straight (sequence of 5)
    if unique_dice == {1,2,3,4,5} or unique_dice == {2,3,4,5,6}:
        scores.append(40)
    else:
        scores.append(0)

    # Yahtzee (five of a kind)
    if 5 in counts.values():
        scores.append(50)
    else:
        scores.append(0)

    # Chance
    scores.append(total_sum)

    return max(scores)

def find_best_move():
    """
    Analyzes the initial hand to find the subset of dice to keep
    that maximizes the expected score after one reroll.
    """
    initial_hand = [3, 3, 3, 5, 6]
    max_expected_score = -1
    best_dice_to_keep = ()
    
    # Use a dictionary to store results for unique sets of kept dice
    # to avoid re-computation (e.g. keeping first two '3's is same as keeping first and third '3').
    analyzed_keeps = {}

    # Iterate through all possible numbers of dice to keep (from 0 to 5)
    for num_to_keep in range(len(initial_hand) + 1):
        # Generate all combinations of dice indices to keep
        for keep_indices in itertools.combinations(range(len(initial_hand)), num_to_keep):
            
            kept_dice_list = [initial_hand[i] for i in keep_indices]
            # Use a sorted tuple as a key for the dictionary
            kept_dice_tuple = tuple(sorted(kept_dice_list))
            
            if kept_dice_tuple in analyzed_keeps:
                continue
            
            num_reroll = 5 - len(kept_dice_tuple)
            
            if num_reroll == 0:
                expected_score = calculate_best_score(kept_dice_list)
            else:
                total_score_for_strategy = 0
                num_outcomes = 0
                # Generate all possible outcomes for the rerolled dice
                for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
                    final_hand = kept_dice_list + list(reroll_outcome)
                    total_score_for_strategy += calculate_best_score(final_hand)
                    num_outcomes += 1
                
                expected_score = total_score_for_strategy / num_outcomes
            
            analyzed_keeps[kept_dice_tuple] = expected_score

            if expected_score > max_expected_score:
                max_expected_score = expected_score
                best_dice_to_keep = kept_dice_tuple

    # Output the result. The problem requests printing each number.
    # We will format it as a comma-separated list of the dice to keep.
    if not best_dice_to_keep:
        print("You should keep no dice and reroll all five.")
    else:
        print(*best_dice_to_keep, sep=", ")

# Run the analysis and print the result
find_best_move()