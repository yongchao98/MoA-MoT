import collections
import itertools

def get_best_score(hand):
    """Calculates the highest possible score for a given 5-dice hand."""
    if not isinstance(hand, list) or len(hand) != 5:
        return 0

    scores = collections.defaultdict(int)
    counts = collections.Counter(hand)
    hand_sum = sum(hand)

    # Upper section
    for i in range(1, 7):
        scores['upper_' + str(i)] = counts[i] * i

    # Lower section
    if 5 in counts.values():
        scores['yahtzee'] = 50
    if 4 in counts.values() or 5 in counts.values():
        scores['four_of_a_kind'] = hand_sum
    if sorted(counts.values()) == [2, 3]:
        scores['full_house'] = 25
    if 3 in counts.values() or 4 in counts.values() or 5 in counts.values():
        scores['three_of_a_kind'] = hand_sum

    unique_dice = set(hand)
    is_large_straight = unique_dice in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]
    is_small_straight = is_large_straight or any(s.issubset(unique_dice) for s in [{1,2,3,4}, {2,3,4,5}, {3,4,5,6}])
    
    if is_large_straight:
        scores['large_straight'] = 40
    if is_small_straight:
        scores['small_straight'] = 30
        
    scores['chance'] = hand_sum

    return max(scores.values()) if scores else 0

def analyze_strategies():
    """Analyzes different strategies and finds the best one."""
    initial_hand = [3, 3, 3, 5, 6]

    strategies = {
        "Keep (3, 3, 3)": ([3, 3, 3], 2),
        "Keep (3, 5, 6)": ([3, 5, 6], 2),
        "Keep (3, 3)": ([3, 3], 3),
        "Keep (5, 6)": ([5, 6], 3),
        "Keep (3, 3, 3, 5, 6)": (initial_hand, 0)
    }

    results = {}
    detailed_breakdown = {}

    for name, (kept_dice, num_reroll) in strategies.items():
        if num_reroll == 0:
            expected_score = get_best_score(kept_dice)
            results[name] = expected_score
            continue

        total_score = 0
        num_outcomes = 6 ** num_reroll
        
        # Store counts of best scores for detailed breakdown
        score_counts = collections.defaultdict(int)
        sum_of_scores = collections.defaultdict(int)
        
        # Iterate over all possible outcomes for the dice being rerolled
        for roll_outcome in itertools.product(range(1, 7), repeat=num_reroll):
            final_hand = kept_dice + list(roll_outcome)
            best_score = get_best_score(final_hand)
            total_score += best_score
            
            score_counts[best_score] += 1
            sum_of_scores[best_score] += best_score

        results[name] = total_score / num_outcomes
        detailed_breakdown[name] = (sum_of_scores, num_outcomes)

    best_strategy_name = max(results, key=results.get)
    
    print("Analysis of reroll strategies for the hand [3, 3, 3, 5, 6]:\n")
    for name, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
        print(f"- {name}: Expected Score = {score:.2f}")

    print("\n" + "="*50)
    print(f"The best strategy is to keep the three 3s.\n")
    print(f"Let's break down the expected score for '{best_strategy_name}':")
    
    sum_of_scores, num_outcomes = detailed_breakdown[best_strategy_name]
    
    equation_parts = []
    total_sum_for_verification = 0
    
    print(f"\nWe reroll 2 dice, leading to {num_outcomes} possible outcomes.")
    print("The final score depends on the outcome of the roll:")

    for score_value, total in sorted(sum_of_scores.items()):
        count = total // score_value
        print(f"- {count:>2} outcomes result in a best score of {score_value}")
        equation_parts.append(f"({count} * {score_value})")
        total_sum_for_verification += total

    final_equation = f"({ ' + '.join(equation_parts) }) / {num_outcomes}"
    
    print("\nThe calculation for the expected score is:")
    print(f"E = {final_equation}")
    
    final_value = total_sum_for_verification / num_outcomes
    print(f"E = {total_sum_for_verification} / {num_outcomes} = {final_value:.4f}")

if __name__ == '__main__':
    analyze_strategies()