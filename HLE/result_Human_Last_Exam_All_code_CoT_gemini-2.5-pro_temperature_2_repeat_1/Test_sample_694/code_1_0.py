import collections
import itertools

def get_best_score(hand):
    """Calculates the maximum possible score for a given 5-dice hand."""
    if not isinstance(hand, list):
        hand = list(hand)
    counts = collections.Counter(hand)
    hand_sum = sum(hand)
    scores = {}

    # Upper section scores
    for i in range(1, 7):
        scores[f'upper_{i}'] = counts.get(i, 0) * i

    # Yahtzee
    scores['yahtzee'] = 50 if 5 in counts.values() else 0

    # Straights
    unique_dice = sorted(list(set(hand)))
    is_ls = False
    if ''.join(map(str, unique_dice)) in ('12345', '23456'):
        is_ls = True
    scores['large_straight'] = 40 if is_ls else 0
    
    is_ss = False
    if not is_ls: # A large straight also qualifies as a small straight
        straights = ["1234", "2345", "3456"]
        hand_str = "".join(map(str,unique_dice))
        for s in straights:
            if s in hand_str:
                is_ss = True
                break
    scores['small_straight'] = 30 if is_ss or is_ls else 0

    # Full House
    scores['full_house'] = 25 if sorted(counts.values()) == [2, 3] else 0
    
    # Of a Kind and Chance
    max_count = 0
    if counts:
      max_count = max(counts.values())
      
    scores['four_of_a_kind'] = hand_sum if max_count >= 4 else 0
    scores['three_of_a_kind'] = hand_sum if max_count >= 3 else 0
    scores['chance'] = hand_sum

    return max(scores.values())

def calculate_expected_score(kept_dice):
    """Calculates the expected score for a strategy of keeping certain dice."""
    num_to_reroll = 5 - len(kept_dice)
    if num_to_reroll == 0:
        return get_best_score(kept_dice)
    
    total_score = 0
    num_outcomes = 6**num_to_reroll
    
    # Iterate through all possible outcomes of the rerolled dice
    for reroll_outcome in itertools.product(range(1, 7), repeat=num_to_reroll):
        final_hand = list(kept_dice) + list(reroll_outcome)
        total_score += get_best_score(final_hand)
        
    return total_score / num_outcomes

def solve():
    """Main function to find the optimal strategy."""
    initial_hand = (3, 3, 3, 5, 6)
    
    # Generate all unique strategies (subsets of dice to keep)
    # Using set to avoid duplicate strategies like keeping (3,5) vs (5,3)
    strategies = set()
    for r in range(len(initial_hand) + 1):
        for combo in itertools.combinations(initial_hand, r):
            strategies.add(tuple(sorted(combo)))
            
    # Calculate expected score for each strategy
    results = {}
    for strategy in sorted(list(strategies), key=len): # Sort for cleaner output
        expected_score = calculate_expected_score(strategy)
        results[strategy] = expected_score
        
    # Find the best strategy
    best_strategy = max(results, key=results.get)
    max_score = results[best_strategy]

    # Print a comparison of the top 5 strategies
    print("--- Yahtzee Strategy Analysis ---")
    print(f"Initial hand: {list(initial_hand)}\n")
    print("Comparison of the top 5 strategies by expected score:")
    
    sorted_results = sorted(results.items(), key=lambda item: item[1], reverse=True)
    
    for strategy, score in sorted_results[:5]:
        strategy_str = f"{list(strategy)}" if strategy else "None"
        print(f"  Keep {strategy_str:<12} | Expected Score = {score:.2f}")

    # Final conclusion based on the calculations
    best_strategy_list = ', '.join(map(str, best_strategy)) if best_strategy else "None"
    print("\n--- Final Recommendation ---")
    print(f"The optimal move is to keep the dice: {best_strategy_list}")
    print(f"This maximizes your expected score for the turn at {max_score:.2f} points.")
    
if __name__ == '__main__':
    solve()
    print("\nTo follow the prompt, the values to keep are formatted below.")
    print("<<<3, 3, 3>>>")