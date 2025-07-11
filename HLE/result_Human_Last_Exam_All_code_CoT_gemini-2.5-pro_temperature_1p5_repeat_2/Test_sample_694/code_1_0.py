import collections
from itertools import product

def calculate_best_score(hand):
    """
    Calculates the highest possible score for a given 5-dice hand.
    """
    counts = collections.Counter(hand)
    total_sum = sum(hand)
    scores = []

    # Fallback score is always the sum of dice (Chance)
    scores.append(total_sum)

    # Lower Section Scores
    # Three/Four of a Kind score is the sum of all dice
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)
    # Yahtzee
    if 5 in counts.values():
        scores.append(50)
    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    
    # Straights
    unique_dice = sorted(list(set(hand)))
    # Large Straight
    if unique_dice in ([1, 2, 3, 4, 5], [2, 3, 4, 5, 6]):
        scores.append(40)
    # Small Straight
    is_small_straight = False
    for i in range(1, 4):
        if set(range(i, i + 4)).issubset(set(unique_dice)):
            is_small_straight = True
            break
    if is_small_straight:
        scores.append(30)
        
    return max(scores)

def calculate_expected_score(dice_to_keep):
    """
    Calculates the expected score for keeping a given set of dice and rerolling the rest.
    """
    num_reroll = 5 - len(dice_to_keep)
    if num_reroll == 0:
        return calculate_best_score(dice_to_keep)

    total_score = 0
    # Generate all possible outcomes for the dice being rerolled
    # The `product` function creates the Cartesian product, like nested loops
    outcomes = product(range(1, 7), repeat=num_reroll)
    num_outcomes = 6 ** num_reroll
    
    for roll in outcomes:
        final_hand = dice_to_keep + list(roll)
        total_score += calculate_best_score(final_hand)
        
    return total_score / num_outcomes

def solve_yahtzee_reroll():
    """
    Analyzes the Yahtzee hand and determines the best dice to keep.
    """
    initial_hand = [3, 3, 3, 5, 6]
    
    # Define the strategies to test
    strategies = {
        "Keep 3, 3, 3": [3, 3, 3],
        "Keep 3, 3, 3, 5": [3, 3, 3, 5],
        "Keep 3, 3, 3, 6": [3, 3, 3, 6],
    }
    
    best_strategy = None
    max_expected_score = -1
    results = {}

    print("Calculating expected scores for different strategies...")
    print("-" * 50)

    for name, kept_dice in strategies.items():
        expected_score = calculate_expected_score(kept_dice)
        results[name] = (kept_dice, expected_score)
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy = kept_dice
            
    # Print the results for each strategy
    for name, (kept_dice, score) in results.items():
        reroll_count = 5 - len(kept_dice)
        print(f"Strategy: {name} (reroll {reroll_count} {'die' if reroll_count == 1 else 'dice'})")
        print(f"Expected Score: {score:.4f}")
        # The following shows the calculation for the expected score
        print(f"Calculation: Sum of best scores for all {6**reroll_count} outcomes / {6**reroll_count}")
        print("-" * 50)
        
    print(f"\nThe best strategy is to keep the dice with values: {', '.join(map(str, best_strategy))}")
    print(f"This yields the highest expected score of {max_expected_score:.4f}")

# Run the analysis
solve_yahtzee_reroll()