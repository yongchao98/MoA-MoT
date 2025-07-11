import collections
from itertools import chain, combinations, product

def get_best_score(dice):
    """
    Calculates the best possible score for a 5-dice hand.
    """
    if not (isinstance(dice, list) or isinstance(dice, tuple)) or len(dice) != 5:
        return 0
    
    counts = collections.Counter(dice)
    total_sum = sum(dice)
    scores = []

    # Upper Section
    for i in range(1, 7):
        scores.append(counts[i] * i)

    # Lower Section
    # Three/Four of a Kind & Chance
    scores.append(total_sum) # Chance is always an option
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum) # 3 of a Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum) # 4 of a Kind

    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)

    # Yahtzee
    if 5 in counts.values():
        scores.append(50)

    # Straights
    unique_dice = sorted(list(set(dice)))
    # Create strings of unique dice to check for sequences
    s = "".join(map(str, unique_dice))
    is_ls = s == "12345" or s == "23456"
    is_ss = "1234" in s or "2345" in s or "3456" in s
    
    if is_ls:
        scores.append(40)
    elif is_ss: # A large straight also counts as a small straight
        scores.append(30)
        
    return max(scores)

def solve_yahtzee_turn():
    """
    Determines the best dice to keep to maximize expected score.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all subsets of indices (0,1,2,3,4) to represent keepers
    indices = list(range(len(initial_hand)))
    powerset = chain.from_iterable(combinations(indices, r) for r in range(len(indices) + 1))
    
    best_strategy = None
    max_expected_score = -1.0
    best_strategy_details = {}

    # Use a set to only check unique sets of keeper values
    checked_strategies = set()

    for keep_indices in powerset:
        dice_to_keep = sorted([initial_hand[i] for i in keep_indices])
        
        strategy_tuple = tuple(dice_to_keep)
        if strategy_tuple in checked_strategies:
            continue
        checked_strategies.add(strategy_tuple)

        num_reroll = 5 - len(dice_to_keep)
        total_score_for_strategy = 0
        
        outcomes = []
        if num_reroll == 0:
            total_score_for_strategy = get_best_score(dice_to_keep)
            outcomes.append({'roll': [], 'final_hand': dice_to_keep, 'score': total_score_for_strategy})
            num_outcomes = 1
        else:
            possible_rolls = product(range(1, 7), repeat=num_reroll)
            num_outcomes = 6 ** num_reroll
            for roll in possible_rolls:
                final_hand = sorted(dice_to_keep + list(roll))
                score = get_best_score(final_hand)
                total_score_for_strategy += score
                outcomes.append({'roll': roll, 'final_hand': final_hand, 'score': score})
        
        expected_score = total_score_for_strategy / num_outcomes
        
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy = dice_to_keep
            best_strategy_details = {
                "total_score": total_score_for_strategy,
                "num_outcomes": num_outcomes,
                "outcomes": outcomes
            }
    
    print(f"The initial hand is: {list(initial_hand)}")
    print("-" * 40)
    print(f"The optimal dice to keep are: {best_strategy}")
    print(f"This yields the highest expected score of: {max_expected_score:.4f}")
    print("-" * 40)
    print(f"This is calculated by considering all {best_strategy_details['num_outcomes']} possible outcomes of rerolling the other dice.")
    
    # Show the "equation" for the best strategy
    score_list = [str(o['score']) for o in best_strategy_details['outcomes']]
    
    # To avoid an unreadably long line, we'll summarize the scores
    score_counts = collections.Counter(score_list)
    equation_parts = [f"({count} outcomes with score {score})" for score, count in sorted(score_counts.items(), key=lambda x: int(x[0]))]

    print("Breakdown of scores for the best strategy:")
    for part in equation_parts:
        print(f"  - {part}")

    total_sum_str = " + ".join([f"{c}*{s}" for s,c in sorted(score_counts.items(), key=lambda x: int(x[0]))])
    
    print(f"\nFinal calculation:")
    print(f"({total_sum_str}) / {best_strategy_details['num_outcomes']} = {best_strategy_details['total_score']} / {best_strategy_details['num_outcomes']} = {max_expected_score:.4f}")


solve_yahtzee_turn()