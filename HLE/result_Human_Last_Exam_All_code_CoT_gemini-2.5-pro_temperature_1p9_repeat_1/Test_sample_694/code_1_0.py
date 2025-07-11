import collections
import itertools

def calculate_best_score(hand):
    """
    Calculates the highest possible score for a given 5-dice hand.
    It checks every Yahtzee category and returns the maximum score.
    """
    counts = collections.Counter(hand)
    total_sum = sum(hand)
    scores = []

    # Upper Section
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Lower Section
    # Three/Four-of-a-kind: Score is sum of all dice
    if max(counts.values()) >= 3:
        scores.append(total_sum)
    if max(counts.values()) >= 4:
        scores.append(total_sum)

    # Full House: 25 points
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)

    # Straights
    unique_dice = set(hand)
    # Small Straight: 30 points
    if {1, 2, 3, 4} <= unique_dice or \
       {2, 3, 4, 5} <= unique_dice or \
       {3, 4, 5, 6} <= unique_dice:
        scores.append(30)

    # Large Straight: 40 points
    if unique_dice == {1, 2, 3, 4, 5} or unique_dice == {2, 3, 4, 5, 6}:
        scores.append(40)

    # Yahtzee: 50 points
    if 5 in counts.values():
        scores.append(50)

    # Chance: Sum of all dice
    scores.append(total_sum)

    return max(scores)

def solve_yahtzee_reroll():
    """
    Calculates the expected score for every possible reroll strategy
    for the given hand and identifies the best option.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    # Use a set to store unique combinations of dice to keep
    unique_keeps = set()
    for r in range(len(initial_hand) + 1):
        for indices in itertools.combinations(range(len(initial_hand)), r):
            # Create a sorted tuple of the dice to keep
            kept_dice = tuple(sorted([initial_hand[i] for i in indices]))
            unique_keeps.add(kept_dice)
    
    strategy_results = {}

    print(f"Analyzing hand: {initial_hand}\n")
    print("Expected score for each possible 'keep' strategy:")
    print("-" * 50)
    # Header for the output table
    print(f"{'Keep Dice':<20} | {'Reroll #':<10} | {'Expected Score':<15}")
    print("-" * 50)

    # Calculate expected score for each unique keep
    for kept_dice in sorted(list(unique_keeps), key=len):
        num_reroll = 5 - len(kept_dice)
        
        if num_reroll == 0:
            # No reroll, score is deterministic
            expected_score = calculate_best_score(list(kept_dice))
        else:
            total_score_from_outcomes = 0
            num_outcomes = 6 ** num_reroll
            
            # Generate all possible outcomes for the dice being rerolled
            for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
                final_hand = list(kept_dice) + list(reroll_outcome)
                total_score_from_outcomes += calculate_best_score(final_hand)
            
            expected_score = total_score_from_outcomes / num_outcomes
            
        strategy_results[kept_dice] = expected_score
        
        # Print each result as it is calculated
        print(f"{str(kept_dice):<20} | {num_reroll:<10} | {expected_score:<15.2f}")

    # Find the best strategy from the results
    best_strategy = max(strategy_results, key=strategy_results.get)
    max_expected_score = strategy_results[best_strategy]

    print("-" * 50)
    print("\nConclusion:")
    print(f"The optimal strategy is to keep the dice: {best_strategy}")
    print(f"This yields the highest expected score of: {max_expected_score:.2f}")

# Run the analysis
solve_yahtzee_reroll()