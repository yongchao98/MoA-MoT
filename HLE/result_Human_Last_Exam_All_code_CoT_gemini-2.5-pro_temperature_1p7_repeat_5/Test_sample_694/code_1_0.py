import itertools
from collections import Counter

def calculate_best_score(dice):
    """Calculates the best possible score for a 5-dice hand."""
    if not isinstance(dice, (list, tuple)) or len(dice) != 5:
        return 0
    
    counts = Counter(dice)
    total_sum = sum(dice)
    scores = []

    # Upper section scores
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Lower section scores
    count_values = list(counts.values())

    # 3-of-a-Kind and 4-of-a-Kind
    scores.append(total_sum if any(c >= 3 for c in count_values) else 0)
    scores.append(total_sum if any(c >= 4 for c in count_values) else 0)

    # Full House
    scores.append(25 if sorted(count_values) == [2, 3] else 0)

    # Small Straight (4 sequential dice)
    unique_dice_set = set(dice)
    s_straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    scores.append(30 if any(s.issubset(unique_dice_set) for s in s_straights) else 0)

    # Large Straight (5 sequential dice)
    scores.append(40 if len(unique_dice_set) == 5 and (max(unique_dice_set) - min(unique_dice_set) == 4) else 0)
    
    # Yahtzee (5 of a kind)
    scores.append(50 if 5 in count_values else 0)

    # Chance
    scores.append(total_sum)
    
    return max(scores)

def solve_yahtzee_reroll():
    """Determines the optimal dice to keep to maximize expected score."""
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate unique subsets of dice to keep
    options_to_keep = set()
    for i in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, i):
            options_to_keep.add(tuple(sorted(combo)))
    
    results = []

    # For each option, calculate the expected score
    for kept_dice in sorted(list(options_to_keep), key=lambda x: (len(x), x)):
        num_reroll = 5 - len(kept_dice)
        
        if num_reroll == 0:
            expected_score = calculate_best_score(kept_dice)
        else:
            total_score = 0
            # Iterate through all possible outcomes of the rerolled dice
            # (e.g., 6^num_reroll outcomes)
            reroll_outcomes = list(itertools.product(range(1, 7), repeat=num_reroll))
            num_outcomes = len(reroll_outcomes)
            
            for reroll_outcome in reroll_outcomes:
                final_hand = kept_dice + reroll_outcome
                total_score += calculate_best_score(final_hand)
                
            expected_score = total_score / num_outcomes
            
        results.append({'kept': kept_dice, 'score': expected_score})

    # Sort results by score to find the best option
    results.sort(key=lambda x: x['score'], reverse=True)
    
    print("Expected score for each possible set of dice to keep:")
    print("-" * 55)
    print(f"{'Dice Kept':<20} | {'Expected Score':<20}")
    print("-" * 55)
    for result in results:
        kept_str = str(result['kept']) if result['kept'] else "None (Reroll all)"
        print(f"{kept_str:<20} | {result['score']:.4f}")

    best_option = results[0]
    best_kept_dice = best_option['kept']
    best_score = best_option['score']
    
    print("\n" + "="*55)
    print("CONCLUSION:")
    print(f"To maximize your expected score, you should keep the dice: {best_kept_dice}")
    print(f"This move yields the highest expected score of {best_score:.4f}.")
    print("="*55)


solve_yahtzee_reroll()