import collections
import itertools

def calculate_best_score(dice):
    """
    Calculates the highest possible score for a given 5-dice hand
    by checking all Yahtzee scoring categories.
    """
    if not isinstance(dice, list):
        dice = list(dice)
    
    counts = collections.Counter(dice)
    total_sum = sum(dice)
    scores = []

    # Upper Section Scores
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Lower Section Scores
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)  # Three of a Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)  # Four of a Kind
    if any(c >= 5 for c in counts.values()):
        scores.append(50)  # Yahtzee
        
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)  # Full House

    # Straights
    unique_dice_set = set(dice)
    # Small Straight (value 30)
    if ({1, 2, 3, 4}.issubset(unique_dice_set) or \
        {2, 3, 4, 5}.issubset(unique_dice_set) or \
        {3, 4, 5, 6}.issubset(unique_dice_set)):
        scores.append(30)
        
    # Large Straight (value 40)
    if len(unique_dice_set) == 5 and (max(dice) - min(dice) == 4):
         scores.append(40)

    scores.append(total_sum)  # Chance

    return max(scores) if scores else 0

def main():
    """
    Analyzes the initial dice to find the optimal keep strategy.
    """
    initial_dice = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    all_possible_keeps = set()
    for i in range(len(initial_dice) + 1):
        for combo in itertools.combinations(initial_dice, i):
            all_possible_keeps.add(tuple(sorted(combo)))

    results = {}

    # Calculate expected score for each subset
    for kept_dice in all_possible_keeps:
        num_to_roll = 5 - len(kept_dice)
        
        if num_to_roll == 0:
            expected_score = calculate_best_score(list(kept_dice))
        else:
            total_score = 0
            num_outcomes = 6**num_to_roll
            # Iterate through all possible reroll outcomes
            for roll_outcome in itertools.product(range(1, 7), repeat=num_to_roll):
                final_dice = list(kept_dice) + list(roll_outcome)
                total_score += calculate_best_score(final_dice)
            expected_score = total_score / num_outcomes
            
        results[kept_dice] = expected_score

    # Find and print the best option
    best_keep = max(results, key=results.get)
    max_score = results[best_keep]
    
    print(f"Initial Dice: {initial_dice}\n")
    print("Based on calculating the expected score for all possibilities, here is the best option:\n")
    print(f"To maximize your expected score, you should keep the dice with the values:")
    
    # Print the numbers in the best set to keep
    keep_str = ", ".join(map(str, best_keep)) if best_keep else "Nothing"
    print(keep_str)
    
    print(f"\nThis strategy results in a maximum expected score of {max_score:.2f} points.")


if __name__ == "__main__":
    main()