import collections

def get_best_score(hand):
    """
    Calculates the best possible score for a given 5-dice hand in Yahtzee.
    Ignores upper section scores as sum-based scores are generally better for this problem.
    """
    scores = [0]
    counts = collections.Counter(hand)
    
    # Chance and N-of-a-Kind (score is sum of all dice)
    total_sum = sum(hand)
    max_count = 0
    if counts:
        max_count = max(counts.values())

    if max_count >= 3:
        scores.append(total_sum) # 3-of-a-kind and 4-of-a-kind score
    else:
        scores.append(total_sum) # Chance score
        
    # Yahtzee
    if max_count == 5:
        scores.append(50)
        
    # Full House
    if 3 in counts.values() and 2 in counts.values():
        scores.append(25)
        
    # Large Straight (sequence of 5)
    unique_dice = sorted(counts.keys())
    if len(unique_dice) == 5 and (unique_dice[-1] - unique_dice[0] == 4):
        scores.append(40)
        
    # Small Straight (sequence of 4)
    s = "".join(map(str, sorted(list(set(hand)))))
    if "1234" in s or "2345" in s or "3456" in s:
        scores.append(30)
        
    return max(scores)

def main():
    initial_hand = [3, 3, 3, 5, 6]
    strategies = {}

    # --- Strategy 1: Keep all dice ---
    score_keep_all = get_best_score(initial_hand)
    strategies['Keep {3, 3, 3, 5, 6}'] = score_keep_all
    
    print("Strategy 1: Keep all dice {3, 3, 3, 5, 6}")
    print(f"The hand {initial_hand} scores as a Three of a Kind.")
    print(f"Score = 3 + 3 + 3 + 5 + 6 = {score_keep_all}\n")

    # --- Strategy 2: Keep {3, 3, 3}, reroll 2 dice ---
    kept_dice = [3, 3, 3]
    total_score = 0
    num_outcomes = 6 * 6
    for d1 in range(1, 7):
        for d2 in range(1, 7):
            final_hand = kept_dice + [d1, d2]
            total_score += get_best_score(final_hand)
    expected_score = total_score / num_outcomes
    strategies['Keep {3, 3, 3}'] = expected_score

    print("Strategy 2: Keep {3, 3, 3}, reroll 2 dice")
    print(f"Sum of best scores over all 36 outcomes = {total_score}")
    print(f"Expected Score = {total_score} / {num_outcomes} = {expected_score:.2f}\n")
    
    # --- Strategy 3: Keep {3, 3, 3, 5}, reroll 1 die ---
    kept_dice = [3, 3, 3, 5]
    total_score = 0
    num_outcomes = 6
    outcomes_str = []
    for d1 in range(1, 7):
        final_hand = kept_dice + [d1]
        score = get_best_score(final_hand)
        total_score += score
        outcomes_str.append(str(score))
    expected_score = total_score / num_outcomes
    strategies['Keep {3, 3, 3, 5}'] = expected_score

    print("Strategy 3: Keep {3, 3, 3, 5}, reroll 1 die")
    print(f"Sum of best scores for each roll = {' + '.join(outcomes_str)} = {total_score}")
    print(f"Expected Score = {total_score} / {num_outcomes} = {expected_score:.2f}\n")

    # --- Strategy 4: Keep {3, 3, 3, 6}, reroll 1 die ---
    kept_dice = [3, 3, 3, 6]
    total_score = 0
    num_outcomes = 6
    outcomes_str = []
    for d1 in range(1, 7):
        final_hand = kept_dice + [d1]
        score = get_best_score(final_hand)
        total_score += score
        outcomes_str.append(str(score))
    expected_score = total_score / num_outcomes
    strategies['Keep {3, 3, 3, 6}'] = expected_score
    
    print("Strategy 4: Keep {3, 3, 3, 6}, reroll 1 die")
    print(f"Sum of best scores for each roll = {' + '.join(outcomes_str)} = {total_score}")
    print(f"Expected Score = {total_score} / {num_outcomes} = {expected_score:.2f}\n")

    # --- Conclusion ---
    best_strategy = max(strategies, key=strategies.get)
    max_score = strategies[best_strategy]
    
    print("--- Conclusion ---")
    print("Comparing the scores from each strategy:")
    for strategy, score in strategies.items():
        print(f"  - {strategy}: {score:.2f}")

    print(f"\nThe best strategy is to '{best_strategy}' with a score of {max_score:.2f}.")
    print("Therefore, you should keep all the dice: 3, 3, 3, 5, 6.")

if __name__ == "__main__":
    main()