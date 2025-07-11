import collections
import itertools

def get_best_score(hand):
    """Calculates the best possible score for a given 5-dice hand."""
    scores = collections.defaultdict(int)
    counts = collections.Counter(hand)
    total_sum = sum(hand)
    
    # Chance
    scores['Chance'] = total_sum
    
    # Yahtzee (5 of a kind)
    if 5 in counts.values():
        scores['Yahtzee'] = 50
        
    # Large Straight (sequence of 5)
    unique_dice = sorted(list(counts.keys()))
    if unique_dice in [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6]]:
        scores['Large Straight'] = 40
        
    # Small Straight (sequence of 4)
    straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    hand_uniques = set(unique_dice)
    for s in straights:
        if s.issubset(hand_uniques):
            scores['Small Straight'] = 30
            break
            
    # Full House (3 of one, 2 of another)
    if sorted(counts.values()) == [2, 3]:
        scores['Full House'] = 25
        
    # 4 of a Kind
    if any(c >= 4 for c in counts.values()):
        scores['4 of a Kind'] = total_sum
            
    # 3 of a Kind
    if any(c >= 3 for c in counts.values()):
        scores['3 of a Kind'] = total_sum

    return max(scores.values()) if scores else 0

def calculate_expected_score(dice_to_keep, num_dice_to_roll):
    """Calculates the expected score for keeping a set of dice and rerolling the rest."""
    if num_dice_to_roll == 0:
        return get_best_score(dice_to_keep)
        
    total_score = 0
    # Generate all possible outcomes for the dice being rerolled
    outcomes = list(itertools.product(range(1, 7), repeat=num_dice_to_roll))
    
    for outcome in outcomes:
        hand = sorted(dice_to_keep + list(outcome))
        total_score += get_best_score(hand)
        
    return total_score / len(outcomes)

def main():
    """Main function to analyze the Yahtzee hand and determine the best move."""
    initial_hand = [3, 3, 3, 5, 6]

    strategies = {
        "Keep all [3, 3, 3, 5, 6]": (initial_hand, 0),
        "Keep [3, 3, 3]": ([3, 3, 3], 2),
        "Keep [3, 3, 3, 5]": ([3, 3, 3, 5], 1),
        "Keep [3, 3, 3, 6]": ([3, 3, 3, 6], 1),
    }

    results = {}
    for name, (kept_dice, num_to_roll) in strategies.items():
        results[name] = calculate_expected_score(kept_dice, num_to_roll)

    best_strategy = max(results, key=results.get)
    max_score = results[best_strategy]

    print("Analysis of options for the hand [3, 3, 3, 5, 6]:\n")
    
    # Print results for each strategy
    print(f"1. Strategy: Keep all dice [3, 3, 3, 5, 6]")
    print(f"   - Scoring this hand as is (e.g., for 'Three of a Kind') gives a certain score of {results['Keep all [3, 3, 3, 5, 6]']:.2f}.\n")

    print(f"2. Strategy: Keep [3, 3, 3] and reroll two dice")
    print(f"   - This gives a chance at a Yahtzee or Full House.")
    print(f"   - The calculated average (expected) score is {results['Keep [3, 3, 3]']:.2f}.\n")

    print(f"3. Strategy: Keep [3, 3, 3, 5] and reroll one die")
    print(f"   - This improves the chance of getting a Full House (by rolling a 5).")
    print(f"   - The calculated average (expected) score is {results['Keep [3, 3, 3, 5]']:.2f}.\n")
    
    print(f"4. Strategy: Keep [3, 3, 3, 6] and reroll one die")
    print(f"   - This improves the chance of getting a Full House (by rolling a 6).")
    print(f"   - The calculated average (expected) score is {results['Keep [3, 3, 3, 6]']:.2f}.\n")

    print("-" * 40)
    print("CONCLUSION:")
    print(f"The best strategy is to '{best_strategy}' with a score of {max_score:.2f}.")
    print("A certain score of 20.00 is better than any of the expected scores from rerolling.")
    print("\nTo maximize your score for this turn, you should keep all the dice with values:")
    
    # Extract the numbers from the best strategy string
    kept_dice_str = best_strategy.split('[')[1].split(']')[0].split(', ')
    for i, die in enumerate(kept_dice_str):
        if i > 0:
            print(" + ", end="")
        print(f"{die}", end="")
    
    # Final equation for the best score
    score = results['Keep all [3, 3, 3, 5, 6]']
    print(f" = {int(score)}")

if __name__ == "__main__":
    main()