import collections
import itertools

def calculate_score(dice):
    """Calculates the best possible score for a given 5-dice hand."""
    counts = collections.Counter(dice)
    total = sum(dice)
    scores = []

    # Lower Section Scoring
    # Yahtzee
    if 5 in counts.values():
        scores.append(50)
    
    # Large Straight
    unique_dice = sorted(counts.keys())
    if unique_dice == [1, 2, 3, 4, 5] or unique_dice == [2, 3, 4, 5, 6]:
        scores.append(40)
    # Small Straight
    else:
        is_ss = False
        # Check for all three possible small straights
        for s in [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]:
            if s.issubset(unique_dice):
                is_ss = True
                break
        if is_ss:
            scores.append(30)
            
    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
        
    # Four of a Kind (score is sum of all dice)
    if 4 in counts.values() or 5 in counts.values():
        scores.append(total)
        
    # Three of a Kind (score is sum of all dice)
    if 3 in counts.values() or 4 in counts.values() or 5 in counts.values():
        scores.append(total)

    # Chance (score is sum of all dice)
    scores.append(total)
    
    return max(scores)

def calculate_expected_score(kept_dice):
    """
    Calculates the expected score for a given set of kept dice by enumerating
    all possible outcomes of a reroll.
    """
    num_reroll = 5 - len(kept_dice)
    if num_reroll == 0:
        return calculate_score(kept_dice)

    total_score = 0
    num_outcomes = 6 ** num_reroll
    
    # Iterate through all possible outcomes of the rerolled dice
    for roll in itertools.product(range(1, 7), repeat=num_reroll):
        final_dice = kept_dice + roll
        total_score += calculate_score(final_dice)
        
    return total_score / num_outcomes

def main():
    """Main function to run the Yahtzee strategy analysis."""
    # Define the strategies (sets of dice to keep) to evaluate
    strategies = {
        "Keep 3, 5, 6 (go for straight)": (3, 5, 6),
        "Keep 3, 3, 5, 6 (pair + straight potential)": (3, 3, 5, 6),
        "Keep 3, 3, 3, 6 (three of a kind + kicker)": (3, 3, 3, 6),
        "Keep 3, 3, 3, 5 (three of a kind + kicker)": (3, 3, 3, 5),
        "Keep 3, 3, 3 (pure three of a kind)": (3, 3, 3),
        "Keep current hand (no reroll)": (3, 3, 3, 5, 6),
        "Keep nothing (reroll all 5 dice)": (),
    }

    # Calculate and store the expected score for each strategy
    results = {}
    print("Analyzing strategies for the hand: 3, 3, 3, 5, 6...")
    for name, kept_dice in strategies.items():
        results[name] = calculate_expected_score(kept_dice)

    # Find the best strategy by the highest expected score
    best_strategy_name = max(results, key=results.get)
    best_score = results[best_strategy_name]
    best_dice_to_keep = strategies[best_strategy_name]

    # Print the analysis results, sorted from best to worst
    print("\n--- Analysis Complete ---")
    print("Expected scores for different strategies (from best to worst):\n")
    for name, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        dice_str = ", ".join(map(str, strategies[name])) if strategies[name] else "None"
        print(f"- To {name}: {score:.2f} points")

    print("\n--- Conclusion ---")
    print(f"The optimal strategy is to '{best_strategy_name}'.")
    final_answer_str = ", ".join(map(str, best_dice_to_keep))
    print(f"You should keep the dice with values: {final_answer_str}")
    print(f"This yields the highest possible expected score of {best_score:.2f} points.")

if __name__ == "__main__":
    main()