import collections
import itertools

def get_best_score(hand):
    """
    Calculates the maximum possible score for a given 5-dice hand
    by evaluating all Yahtzee scoring categories.
    """
    counts = collections.Counter(hand)
    
    # Upper Section Scores
    scores = [counts[i] * i for i in range(1, 7)]
    
    # Lower Section Scores
    # Three and Four of a Kind
    if any(c >= 3 for c in counts.values()):
        scores.append(sum(hand))
    if any(c >= 4 for c in counts.values()):
        scores.append(sum(hand))
    
    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    
    # Straights
    unique_dice = set(hand)
    is_small_straight = (
        {1, 2, 3, 4}.issubset(unique_dice) or
        {2, 3, 4, 5}.issubset(unique_dice) or
        {3, 4, 5, 6}.issubset(unique_dice)
    )
    is_large_straight = (
        unique_dice == {1, 2, 3, 4, 5} or
        unique_dice == {2, 3, 4, 5, 6}
    )
    if is_small_straight:
        scores.append(30)
    if is_large_straight:
        scores.append(40)
    
    # Yahtzee
    if 5 in counts.values():
        scores.append(50)
    
    # Chance
    scores.append(sum(hand))
    
    return max(scores) if scores else 0

def find_best_strategy():
    """
    Analyzes all possible strategies for the given hand and identifies the optimal one.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    # We use indices to correctly handle duplicate dice values
    hand_indices = list(range(len(initial_hand)))
    strategies = set()
    for i in range(len(hand_indices) + 1):
        for subset_indices in itertools.combinations(hand_indices, i):
            dice_to_keep = tuple(sorted([initial_hand[j] for j in subset_indices]))
            strategies.add(dice_to_keep)

    best_strategy = None
    max_expected_score = -1

    for kept_dice in strategies:
        num_to_reroll = 5 - len(kept_dice)
        
        if num_to_reroll == 0:
            expected_score = get_best_score(kept_dice)
        else:
            total_score_for_strategy = 0
            num_outcomes = 6 ** num_to_reroll
            
            # Iterate through all possible outcomes of the reroll
            for roll_outcome in itertools.product(range(1, 7), repeat=num_to_reroll):
                final_hand = kept_dice + roll_outcome
                total_score_for_strategy += get_best_score(final_hand)
            
            expected_score = total_score_for_strategy / num_outcomes
            
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy = kept_dice

    # The problem asks to output the values to keep
    final_values = ', '.join(map(str, best_strategy))
    if not final_values:
        final_values = "None (reroll all dice)"

    print("After analyzing all possible sets of dice to keep, the optimal strategy is:")
    print(f"Values to keep: {final_values}")
    print(f"This strategy yields the highest expected score of approximately {max_expected_score:.2f} points.")
    
# Run the analysis and print the result
find_best_strategy()