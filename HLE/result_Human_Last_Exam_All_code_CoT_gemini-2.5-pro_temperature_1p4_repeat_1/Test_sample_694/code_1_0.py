import collections
import itertools
from math import floor

def get_hand_counts(hand):
    """Uses collections.Counter to get the counts of each die face."""
    return collections.Counter(hand)

# --- Yahtzee Scoring Functions ---

def score_upper_section(hand, num):
    """Calculates the score for a specific number in the upper section."""
    return hand.count(num) * num

def score_n_of_a_kind(hand, n):
    """Calculates score for Three and Four of a Kind."""
    counts = get_hand_counts(hand)
    if any(c >= n for c in counts.values()):
        return sum(hand)
    return 0

def score_full_house(hand):
    """Calculates score for a Full House."""
    counts = get_hand_counts(hand)
    if sorted(counts.values()) == [2, 3]:
        return 25
    return 0

def score_small_straight(hand):
    """Calculates score for a Small Straight."""
    # Normalize to a set of unique numbers
    hand_set = set(hand)
    # Check for all possible small straights
    if {1, 2, 3, 4}.issubset(hand_set) or \
       {2, 3, 4, 5}.issubset(hand_set) or \
       {3, 4, 5, 6}.issubset(hand_set):
        return 30
    return 0

def score_large_straight(hand):
    """Calculates score for a Large Straight."""
    hand_set = set(hand)
    if hand_set == {1, 2, 3, 4, 5} or hand_set == {2, 3, 4, 5, 6}:
        return 40
    return 0

def score_yahtzee(hand):
    """Calculates score for a Yahtzee."""
    counts = get_hand_counts(hand)
    if 5 in counts.values():
        return 50
    return 0

def score_chance(hand):
    """Calculates score for Chance."""
    return sum(hand)

def calculate_max_score(hand):
    """Calculates the maximum possible score for a given 5-dice hand."""
    scores = [
        score_upper_section(hand, 1),
        score_upper_section(hand, 2),
        score_upper_section(hand, 3),
        score_upper_section(hand, 4),
        score_upper_section(hand, 5),
        score_upper_section(hand, 6),
        score_n_of_a_kind(hand, 3),
        score_n_of_a_kind(hand, 4),
        score_full_house(hand),
        score_small_straight(hand),
        score_large_straight(hand),
        score_yahtzee(hand),
        score_chance(hand)
    ]
    return max(scores)

def calculate_expected_score(kept_dice):
    """Calculates the expected score for a given set of kept dice."""
    num_reroll = 5 - len(kept_dice)
    
    if num_reroll == 0:
        return calculate_max_score(kept_dice)
    
    total_score = 0
    num_outcomes = 6 ** num_reroll
    
    # Iterate through all possible outcomes of the reroll
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
        final_hand = kept_dice + list(reroll_outcome)
        total_score += calculate_max_score(final_hand)
        
    return total_score / num_outcomes

def solve_yahtzee_reroll():
    """Main function to find the best strategy."""
    initial_hand = [3, 3, 3, 5, 6]
    
    # Generate all unique strategies for keeping dice
    strategies = set()
    for i in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, i):
            strategies.add(tuple(sorted(combo)))

    strategy_scores = {}
    print("Calculating expected score for each possible strategy...")
    print("-" * 50)
    
    for strategy in sorted(list(strategies), key=len):
        score = calculate_expected_score(list(strategy))
        strategy_scores[strategy] = score
        # The prompt asks to "output each number in the final equation"
        # We interpret this as showing the results for each strategy.
        strategy_str = f"Keeping {str(strategy):<20}"
        if not strategy:
            strategy_str = f"Keeping {str('()'):<20}" # Handle empty tuple case
        
        print(f"{strategy_str} -> Expected Score: {score:.4f}")

    # Find the best strategy
    best_strategy = max(strategy_scores, key=strategy_scores.get)
    max_score = strategy_scores[best_strategy]
    
    print("-" * 50)
    print("Conclusion:")
    print(f"The best strategy is to keep the dice: {best_strategy}")
    print(f"This strategy yields the highest expected score of: {max_score:.4f}")
    
    # Format final answer as requested
    final_answer_str = ", ".join(map(str, best_strategy)) if best_strategy else ""
    print("\nBased on this calculation, the values you should keep are:")
    print(final_answer_str)
    # The final answer block is below for programmatic parsing.
    # The printed output above provides the human-readable explanation.

solve_yahtzee_reroll()
<<<3, 3, 3, 5, 6>>>