import itertools
from collections import Counter

def solve_yahtzee_reroll():
    """
    Analyzes a Yahtzee hand to find the optimal dice to keep to maximize
    the expected score after one reroll.
    """
    
    # Memoization cache to store results for already computed scenarios
    memo_expected_score = {}

    def get_score(dice):
        """Calculates the best possible Yahtzee score for a given 5-dice hand."""
        counts = Counter(dice)
        total_sum = sum(dice)
        scores = []

        # Upper section scores
        for i in range(1, 7):
            scores.append(i * counts[i])

        # Chance score
        scores.append(total_sum)

        # Lower section scores
        if any(c >= 3 for c in counts.values()):
            scores.append(total_sum)  # Three of a Kind
        if any(c >= 4 for c in counts.values()):
            scores.append(total_sum)  # Four of a Kind
        if 3 in counts.values() and 2 in counts.values():
            scores.append(25)  # Full House
        if 5 in counts.values():
            scores.append(50)  # Yahtzee

        unique_dice_keys = sorted(counts.keys())
        s = "".join(map(str, unique_dice_keys))
        if "1234" in s or "2345" in s or "3456" in s:
            scores.append(30)  # Small Straight
        if unique_dice_keys == [1, 2, 3, 4, 5] or unique_dice_keys == [2, 3, 4, 5, 6]:
            scores.append(40)  # Large Straight
            
        return max(scores) if scores else 0

    def get_expected_score(kept_dice):
        """
        Calculates the expected best score for a given set of kept dice by
        simulating all possible reroll outcomes.
        """
        kept_dice_key = tuple(sorted(kept_dice))
        if kept_dice_key in memo_expected_score:
            return memo_expected_score[kept_dice_key]

        num_reroll = 5 - len(kept_dice)
        if num_reroll == 0:
            return get_score(kept_dice)

        total_score = 0
        num_outcomes = 6**num_reroll
        
        # Generate all possible outcomes for the dice to be rerolled
        possible_rolls = itertools.product(range(1, 7), repeat=num_reroll)

        for roll in possible_rolls:
            final_hand = list(kept_dice) + list(roll)
            total_score += get_score(final_hand)
        
        expected_score = total_score / num_outcomes
        memo_expected_score[kept_dice_key] = expected_score
        return expected_score

    # The user's hand after the first reroll
    initial_hand = [3, 3, 3, 5, 6]

    # Generate all unique subsets of the initial hand to evaluate
    hand_indices = list(range(len(initial_hand)))
    kept_dice_options = set()
    for i in range(len(initial_hand) + 1):
        for indices in itertools.combinations(hand_indices, i):
            subset = tuple(sorted([initial_hand[j] for j in indices]))
            kept_dice_options.add(subset)

    max_expected_score = -1
    best_dice_to_keep = []

    # Calculate expected score for each possibility and find the best one
    for kept in kept_dice_options:
        kept_list = list(kept)
        expected_score = get_expected_score(kept_list)
        
        if expected_score > max_expected_score:
            max_expected_score = expected_score
            best_dice_to_keep = sorted(kept_list)

    print("To maximize the expected score, you should keep the dice with the following values:")
    # The final answer is the set of numbers to keep.
    # The following line prints each number separated by a space.
    print(*best_dice_to_keep)

solve_yahtzee_reroll()