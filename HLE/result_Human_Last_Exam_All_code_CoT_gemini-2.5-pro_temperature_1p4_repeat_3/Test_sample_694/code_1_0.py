import itertools
from collections import Counter

def calculate_best_score(hand):
    """
    Calculates the best possible score for a given 5-dice hand in Yahtzee.
    """
    scores = []
    total_sum = sum(hand)
    
    # 1. Chance is always an option
    scores.append(total_sum)

    # 2. Check categories based on dice counts
    counts = Counter(hand)
    vals = list(counts.values())
    
    # Yahtzee (50 points)
    if 5 in vals:
        scores.append(50)
    
    # Four of a Kind (sum of all dice)
    if 4 in vals or 5 in vals:
        scores.append(total_sum)
        
    # Three of a Kind (sum of all dice)
    if 3 in vals or 4 in vals or 5 in vals:
        scores.append(total_sum)

    # Full House (25 points)
    if sorted(vals) == [2, 3]:
        scores.append(25)

    # 3. Check for Straights
    unique_dice = sorted(list(set(hand)))
    
    # Large Straight (40 points)
    if unique_dice in ([1, 2, 3, 4, 5], [2, 3, 4, 5, 6]):
        scores.append(40)
    
    # Small Straight (30 points)
    is_small_straight = False
    # Check for any sequence of 4 consecutive dice
    for i in range(len(unique_dice) - 3):
        if unique_dice[i+1] == unique_dice[i] + 1 and \
           unique_dice[i+2] == unique_dice[i] + 2 and \
           unique_dice[i+3] == unique_dice[i] + 3:
            is_small_straight = True
            break
    if is_small_straight:
        scores.append(30)
        
    return max(scores)

def get_score_category(hand, score):
    """
    Helper function to determine the name of the category for a given score.
    This is for descriptive output only.
    """
    if score == 50: return "Yahtzee"
    if score == 40: return "Large Straight"
    if score == 30: return "Small Straight"
    if score == 25: return "Full House"
    counts = Counter(hand)
    vals = list(counts.values())
    if 4 in vals or 5 in vals: return "4 of a Kind"
    if 3 in vals: return "3 of a Kind"
    return "Chance"


def analyze_strategies():
    """
    Analyzes different strategies for the Yahtzee hand [3, 3, 3, 5, 6].
    """
    initial_hand = [3, 3, 3, 5, 6]
    
    strategies = [
        ("Keep {3, 3, 3}", [3, 3, 3]),
        ("Keep {3, 3, 3, 5}", [3, 3, 3, 5]),
        ("Keep {3, 3, 3, 6}", [3, 3, 3, 6]),
    ]
    
    results = []

    print(f"Initial hand: {initial_hand}\n")

    for name, dice_to_keep in strategies:
        num_reroll = 5 - len(dice_to_keep)
        num_outcomes = 6 ** num_reroll
        total_score = 0
        
        print(f"--- Analyzing Strategy: {name} ---")
        
        # Generate all possible outcomes for the dice being rerolled
        outcomes = itertools.product(range(1, 7), repeat=num_reroll)
        
        score_list = []
        for outcome in outcomes:
            final_hand = dice_to_keep + list(outcome)
            final_hand.sort()
            score = calculate_best_score(final_hand)
            score_list.append(score)
            total_score += score
            
            if num_reroll == 1:
                category = get_score_category(final_hand, score)
                print(f"  Rerolling a {outcome[0]} -> Final Hand: {final_hand}, Best Score: {score} ({category})")

        expected_value = total_score / num_outcomes
        results.append({'strategy': name, 'dice_kept': dice_to_keep, 'expected_score': expected_value, 'scores': score_list})
        
        # Build and print the equation for the expected value
        equation_str = " + ".join(map(str, score_list))
        print(f"\nExpected Score = ({equation_str}) / {num_outcomes} = {expected_value:.2f}\n")

    # Find and announce the best strategy
    best_strategy = max(results, key=lambda x: x['expected_score'])
    
    print("--- Comparison of Strategies ---")
    for res in results:
        print(f"Strategy: {res['strategy']:<20} | Expected Score: {res['expected_score']:.4f}")
    
    print("\n--- Conclusion ---")
    print(f"The best strategy is to keep the dice {{{', '.join(map(str, sorted(best_strategy['dice_kept'])))}}} and reroll the others.")
    print(f"This strategy yields the highest expected score of {best_strategy['expected_score']:.4f}.")
    
    # Return the answer in the required format
    best_dice_str = ', '.join(map(str, sorted(best_strategy['dice_kept'])))
    return f"<<<{best_dice_str}>>>"

# Run the analysis and print the final answer
final_answer = analyze_strategies()
print(final_answer)
