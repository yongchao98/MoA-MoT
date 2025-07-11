import itertools
from collections import Counter

def get_best_score(hand):
    """Calculates the best possible score for a given 5-dice hand."""
    if not hand:
        return 0
    
    counts = Counter(hand)
    total_sum = sum(hand)
    scores = []

    # Upper Section (e.g., score of all 3s)
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # Lower Section
    # Yahtzee (5 of a kind)
    if 5 in counts.values():
        scores.append(50)

    # Large Straight (sequence of 5)
    unique_dice_sorted = sorted(list(counts.keys()))
    if unique_dice_sorted in [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6]]:
        scores.append(40)

    # Small Straight (sequence of 4)
    unique_dice_set = set(counts.keys())
    if {1, 2, 3, 4}.issubset(unique_dice_set) or \
       {2, 3, 4, 5}.issubset(unique_dice_set) or \
       {3, 4, 5, 6}.issubset(unique_dice_set):
        scores.append(30)

    # Full House (3 of one, 2 of another)
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)

    # Four of a Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)

    # Three of a Kind
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)

    # Chance (sum of all dice)
    scores.append(total_sum)

    return max(scores)

def solve_yahtzee_reroll():
    """
    Determines the optimal dice to keep to maximize expected score in Yahtzee.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    unique_keeps = set()
    for i in range(len(initial_hand) + 1):
        for subset_indices in itertools.combinations(range(len(initial_hand)), i):
            # Sort to treat combinations of the same numbers as identical
            kept_dice = tuple(sorted([initial_hand[j] for j in subset_indices]))
            unique_keeps.add(kept_dice)

    max_expected_score = -1
    best_keep_decision = None

    # Calculate expected score for each possible keep decision
    for kept_dice in sorted(list(unique_keeps), key=len):
        num_reroll = 5 - len(kept_dice)
        
        if num_reroll == 0:
            expected_score = get_best_score(kept_dice)
        else:
            total_score = 0
            num_outcomes = 6 ** num_reroll
            
            # Iterate through all possible outcomes of the reroll
            for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
                final_hand = kept_dice + reroll_outcome
                total_score += get_best_score(final_hand)
            
            expected_score = total_score / num_outcomes

        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_keep_decision = kept_dice

    print("To maximize the expected score, you should keep the dice with the following values:")
    if not best_keep_decision:
        print("None (reroll all dice)")
    else:
        # Print each number of the dice combination to keep
        print(*best_keep_decision)

solve_yahtzee_reroll()
<<<3, 3, 3, 5, 6>>>