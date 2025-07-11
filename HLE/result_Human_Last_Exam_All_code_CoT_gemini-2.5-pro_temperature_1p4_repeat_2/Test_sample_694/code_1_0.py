import itertools
from collections import Counter

def get_best_score(hand):
    """Calculates the maximum possible score for a given 5-dice hand."""
    scores = []
    hand_counts = Counter(hand)
    hand_sum = sum(hand)
    
    # Add scores for upper section (Aces through Sixes)
    for i in range(1, 7):
        scores.append(hand_counts.get(i, 0) * i)
        
    # Add score for Chance
    scores.append(hand_sum)
    
    # Check for lower section categories
    counts = sorted(hand_counts.values(), reverse=True)
    
    # Three-Of-A-Kind (if at least 3 same dice)
    if counts and counts[0] >= 3:
        scores.append(hand_sum)
        
    # Four-Of-A-Kind (if at least 4 same dice)
    if counts and counts[0] >= 4:
        scores.append(hand_sum)
        
    # Yahtzee (Five-Of-A-Kind)
    if counts and counts[0] >= 5:
        scores.append(50)
        
    # Full House (3 of one kind, 2 of another)
    if counts == [3, 2]:
        scores.append(25)
        
    # Large Straight (sequence of 5)
    hand_set = set(hand)
    if hand_set in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]:
        scores.append(40)
        
    # Small Straight (sequence of 4)
    straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    for s in straights:
        if s.issubset(hand_set):
            scores.append(30)
            break
            
    return max(scores)

def calculate_expected_value(kept_dice):
    """Calculates the expected score for a given set of kept dice."""
    num_reroll = 5 - len(kept_dice)
    
    # If no dice are rerolled, the score is fixed.
    if num_reroll == 0:
        return get_best_score(list(kept_dice))
        
    total_score_sum = 0
    # Generate all possible outcomes for the dice being rerolled.
    possible_rolls = list(itertools.product(range(1, 7), repeat=num_reroll))
    num_outcomes = len(possible_rolls)
    
    for roll in possible_rolls:
        final_hand = list(kept_dice) + list(roll)
        total_score_sum += get_best_score(final_hand)
        
    return total_score_sum / num_outcomes

def solve_yahtzee_hand():
    """
    Analyzes the Yahtzee hand and determines the optimal strategy.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    print(f"Analyzing Yahtzee hand: {initial_hand}")
    print("Finding the combination of dice to keep that maximizes the expected score...")
    
    # Generate all unique combinations of dice to keep
    keep_options = set()
    for i in range(len(initial_hand) + 1):
        # Generate combinations of dice from the initial hand
        for combo in itertools.combinations(initial_hand, i):
            # Sort to treat e.g. (3,5) and (5,3) as the same choice
            keep_options.add(tuple(sorted(combo)))
            
    # Calculate and store the expected value for each option
    results = {}
    print("\n--- Analysis of Possible Choices ---")
    # Sort options by the number of dice kept for readability
    sorted_options = sorted(list(keep_options), key=len)

    for option in sorted_options:
        ev = calculate_expected_value(option)
        results[option] = ev
        # Use a descriptive label for keeping nothing
        option_str = str(option) if option else "Nothing (Reroll all)"
        print(f"Keeping {option_str:<22} -> Expected Score: {ev:.4f}")

    # Find the best option by finding the max expected value
    best_option = max(results, key=results.get)
    max_ev = results[best_option]
    
    print("\n--- Conclusion ---")
    print("To maximize your score this turn, you should choose the option with the highest expected value.")
    print(f"The best combination to keep is {best_option} with an expected score of {max_ev:.4f}.")
    print("\nThis is the optimal move because keeping the three 3s gives you a great foundation for several high-scoring categories:")
    print("- A good chance at a Full House (25 points).")
    print("- A possibility of a Four-of-a-Kind or even a Yahtzee (50 points).")
    print("- A guaranteed Three-of-a-Kind for a solid fallback score.")

    print("\nFinal equation for the best choice:")
    equation_str = " + ".join([f"Keep {d}" for d in best_option])
    print(f"{equation_str} + Reroll {5 - len(best_option)} dice = Expected Score of {max_ev:.4f}")

if __name__ == '__main__':
    solve_yahtzee_hand()