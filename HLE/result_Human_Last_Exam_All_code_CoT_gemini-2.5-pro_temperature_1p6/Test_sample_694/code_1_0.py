import collections
import itertools

def calculate_scores(dice):
    """Calculates the score for a given hand in all Yahtzee categories."""
    scores = {}
    counts = collections.Counter(dice)
    total_sum = sum(dice)

    # Upper section
    for i in range(1, 7):
        scores[f'Upper {i}s'] = counts.get(i, 0) * i

    # Lower section
    scores['Chance'] = total_sum

    # N-of-a-kind
    num_most_common = 0
    if counts:
        num_most_common = counts.most_common(1)[0][1]

    if num_most_common >= 3:
        scores['Three-of-a-kind'] = total_sum
    else:
        scores['Three-of-a-kind'] = 0

    if num_most_common >= 4:
        scores['Four-of-a-kind'] = total_sum
    else:
        scores['Four-of-a-kind'] = 0

    if num_most_common == 5:
        scores['Yahtzee'] = 50
    else:
        scores['Yahtzee'] = 0

    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores['Full House'] = 25
    else:
        scores['Full House'] = 0

    # Straights
    unique_dice = set(dice)
    is_ls = (unique_dice == {1, 2, 3, 4, 5}) or (unique_dice == {2, 3, 4, 5, 6})
    is_ss = is_ls or ({1, 2, 3, 4} <= unique_dice) or \
            ({2, 3, 4, 5} <= unique_dice) or ({3, 4, 5, 6} <= unique_dice)
    
    scores['Large Straight'] = 40 if is_ls else 0
    scores['Small Straight'] = 30 if is_ss else 0
    
    return scores

def evaluate_strategy(kept_dice):
    """Calculates the maximum expected score for keeping a set of dice."""
    num_to_roll = 5 - len(kept_dice)
    if num_to_roll == 0:
        return max(calculate_scores(kept_dice).values())

    total_max_score = 0
    num_outcomes = 6 ** num_to_roll

    # Use itertools.product to get all possible dice rolls
    for roll_outcome in itertools.product(range(1, 7), repeat=num_to_roll):
        final_hand = kept_dice + roll_outcome
        scores = calculate_scores(final_hand)
        max_score = max(scores.values())
        total_max_score += max_score
        
    return total_max_score / num_outcomes

def main():
    """Main function to solve the yahtzee problem."""
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    strategies_to_evaluate = set()
    for i in range(len(initial_hand) + 1):
        for subset in itertools.combinations(initial_hand, i):
            strategies_to_evaluate.add(tuple(sorted(subset)))

    print(f"Analyzing hand: {initial_hand}\n")
    results = {}
    for kept_dice in sorted(list(strategies_to_evaluate), key=lambda x: (len(x), x), reverse=True):
        expected_score = evaluate_strategy(kept_dice)
        results[kept_dice] = expected_score
    
    # Find the best strategy
    best_strategy_dice = max(results, key=results.get)
    max_expected_score = results[best_strategy_dice]

    print("--- Results ---")
    for kept, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        dice_str = str(kept) if kept else "nothing"
        print(f"Keeping {dice_str:<15}: Expected Score = {score:.2f}")

    print("\n--- Recommendation ---")
    print(f"To maximize your score this turn, you should keep the dice: {', '.join(map(str, best_strategy_dice))}")
    print(f"This strategy yields the highest expected score of {max_expected_score:.2f} points.")
    
    final_answer = ', '.join(map(str, best_strategy_dice))
    # This is where the final answer for the specialized format is generated
    # The printed output is for the user's benefit
    # The actual final answer is just the dice to keep
    # print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    main()