import collections
import itertools

def get_best_score(dice):
    """Calculates the best possible score for a given hand of 5 dice."""
    if not isinstance(dice, collections.Counter):
        counts = collections.Counter(dice)
    else:
        counts = dice
    
    scores = []
    
    # Upper Section
    for i in range(1, 7):
        scores.append(counts[i] * i)
        
    # Lower Section
    total_sum = sum(dice)
    scores.append(total_sum) # Chance

    # N-of-a-kind
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum) # 3-of-a-Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum) # 4-of-a-Kind
    
    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)
        
    # Straights
    unique_dice = set(dice)
    is_small_straight = (
        {1, 2, 3, 4}.issubset(unique_dice) or
        {2, 3, 4, 5}.issubset(unique_dice) or
        {3, 4, 5, 6}.issubset(unique_dice)
    )
    is_large_straight = (
        unique_dice == {1, 2, 3, 4, 5} or
        unique_dice == {2, 3, 4, 5, 6}
    )
    
    if is_large_straight:
        scores.append(40)
    elif is_small_straight:
        scores.append(30)
        
    # Yahtzee
    if 5 in counts.values():
        scores.append(50)
        
    return max(scores)

def main():
    """
    Analyzes the Yahtzee hand {3, 3, 3, 5, 6} to find the optimal dice to keep.
    """
    initial_dice = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    options_to_keep = set()
    for r in range(len(initial_dice) + 1):
        for indices in itertools.combinations(range(len(initial_dice)), r):
            subset = tuple(sorted([initial_dice[i] for i in indices]))
            options_to_keep.add(subset)

    max_expected_score = 0
    best_option_to_keep = None

    for keep in sorted(list(options_to_keep)):
        num_to_reroll = 5 - len(keep)
        
        if num_to_reroll == 0:
            # No dice to reroll, score is fixed
            expected_score = get_best_score(list(keep))
        else:
            total_score = 0
            num_outcomes = 6 ** num_to_reroll
            
            # Iterate through all possible outcomes of the reroll
            for reroll in itertools.product(range(1, 7), repeat=num_to_reroll):
                final_hand = list(keep) + list(reroll)
                total_score += get_best_score(final_hand)
            
            expected_score = total_score / num_outcomes
            
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_option_to_keep = keep

    print(f"The initial hand is: {', '.join(map(str, initial_dice))}")
    print("\nAnalysis of options:")
    
    # Print a few key comparisons
    keep_all_score = get_best_score(list(initial_dice))
    print(f"Keeping all dice ({', '.join(map(str, initial_dice))}) gives a guaranteed score of: {keep_all_score:.2f}")

    keep_threes_reroll_two = [keep for keep in options_to_keep if keep == (3, 3, 3)][0]
    num_to_reroll = 5 - len(keep_threes_reroll_two)
    total_score = sum(get_best_score(list(keep_threes_reroll_two) + list(reroll)) for reroll in itertools.product(range(1, 7), repeat=num_to_reroll))
    ev_threes = total_score / (6 ** num_to_reroll)
    print(f"Keeping ({', '.join(map(str, keep_threes_reroll_two))}) and rerolling 2 dice gives an expected score of: {ev_threes:.2f}")

    keep_straight_prospects = [keep for keep in options_to_keep if keep == (3, 5, 6)][0]
    num_to_reroll = 5 - len(keep_straight_prospects)
    total_score = sum(get_best_score(list(keep_straight_prospects) + list(reroll)) for reroll in itertools.product(range(1, 7), repeat=num_to_reroll))
    ev_straight = total_score / (6 ** num_to_reroll)
    print(f"Keeping ({', '.join(map(str, keep_straight_prospects))}) and rerolling 2 dice gives an expected score of: {ev_straight:.2f}")

    print("\n---Conclusion---")
    print(f"The best option is to keep the dice: {', '.join(map(str, best_option_to_keep))}")
    print(f"This maximizes the expected score to: {max_expected_score:.2f}")

if __name__ == '__main__':
    main()