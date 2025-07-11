import itertools
from collections import Counter

def score_upper_section(hand, number):
    """Calculates the score for the upper section (Ones to Sixes)."""
    return sum(d for d in hand if d == number)

def score_three_of_a_kind(hand):
    """Calculates the score for Three-of-a-Kind."""
    counts = Counter(hand)
    if any(c >= 3 for c in counts.values()):
        return sum(hand)
    return 0

def score_four_of_a_kind(hand):
    """Calculates the score for Four-of-a-Kind."""
    counts = Counter(hand)
    if any(c >= 4 for c in counts.values()):
        return sum(hand)
    return 0

def score_full_house(hand):
    """Calculates the score for a Full House."""
    counts = Counter(hand)
    if sorted(counts.values()) == [2, 3]:
        return 25
    return 0

def score_small_straight(hand):
    """Calculates the score for a Small Straight."""
    unique_dice = set(hand)
    # Check for all possible 4-dice sequences
    if {1, 2, 3, 4}.issubset(unique_dice) or \
       {2, 3, 4, 5}.issubset(unique_dice) or \
       {3, 4, 5, 6}.issubset(unique_dice):
        return 30
    return 0

def score_large_straight(hand):
    """Calculates the score for a Large Straight."""
    unique_dice = set(hand)
    if unique_dice == {1, 2, 3, 4, 5} or unique_dice == {2, 3, 4, 5, 6}:
        return 40
    return 0

def score_yahtzee(hand):
    """Calculates the score for a Yahtzee."""
    if len(set(hand)) == 1:
        return 50
    return 0

def score_chance(hand):
    """Calculates the score for Chance."""
    return sum(hand)

def get_best_score(hand):
    """Finds the maximum possible score for a given hand."""
    scores = [
        score_upper_section(hand, 1),
        score_upper_section(hand, 2),
        score_upper_section(hand, 3),
        score_upper_section(hand, 4),
        score_upper_section(hand, 5),
        score_upper_section(hand, 6),
        score_three_of_a_kind(hand),
        score_four_of_a_kind(hand),
        score_full_house(hand),
        score_small_straight(hand),
        score_large_straight(hand),
        score_yahtzee(hand),
        score_chance(hand),
    ]
    return max(scores)

def solve_yahtzee_reroll():
    """Main function to analyze the Yahtzee hand and find the best move."""
    initial_hand = [3, 3, 3, 5, 6]
    
    # Generate all unique combinations of dice to keep
    options_to_keep = set()
    for i in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, i):
            options_to_keep.add(tuple(sorted(combo)))
            
    results = {}
    
    # Calculate expected score for each option
    for kept_dice_tuple in sorted(list(options_to_keep), key=len):
        kept_dice = list(kept_dice_tuple)
        num_to_reroll = 5 - len(kept_dice)
        
        if num_to_reroll == 0:
            expected_score = float(get_best_score(kept_dice))
        else:
            possible_rolls = itertools.product(range(1, 7), repeat=num_to_reroll)
            total_outcomes = 6 ** num_to_reroll
            total_score = 0
            for roll in possible_rolls:
                final_hand = kept_dice + list(roll)
                total_score += get_best_score(final_hand)
            expected_score = total_score / total_outcomes
            
        results[kept_dice_tuple] = expected_score
        
    # Print results
    print("Expected score for each possible set of dice to keep:")
    for option, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        option_str = str(option) if option else "Nothing"
        print(f"- Keep {option_str}: {score:.4f} average score")
        
    # Find and print the best option
    best_option = max(results, key=results.get)
    best_score = results[best_option]
    
    print("\n--- Conclusion ---")
    print(f"The best option is to keep the dice: {best_option}")
    print(f"This yields the highest expected score of {best_score:.4f}.")
    
    # Required for the final answer format
    global final_answer
    final_answer = ", ".join(map(str, best_option))


# This is a global variable to store the final answer string for the '<<<...>>>' block
final_answer = ""
solve_yahtzee_reroll()
print(f"\n<<<{final_answer}>>>")