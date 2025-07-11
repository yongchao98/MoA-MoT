import collections

memo = {}

def can_form_sum(target, scores, k):
    """
    Checks if a target sum can be formed by choosing k numbers from the scores list.
    Uses recursion with memoization (dynamic programming) for efficiency.
    """
    state = (target, k)
    if state in memo:
        return memo[state]
    if target == 0 and k == 0:
        return True
    if target < 0 or k <= 0:
        return False
    
    # scores is a tuple, so it's hashable and can be part of the memo key if needed,
    # but it's constant for each sub-problem tree, so it's omitted from the key.
    for s in scores:
        if can_form_sum(target - s, scores, k - 1):
            memo[state] = True
            return True
            
    memo[state] = False
    return False

def solve():
    """
    Main function to find the number of possible bull's eye scores.
    """
    # Let v_i = score_i / 5.
    # Anna's equation: 3*v1 + v2 + v_bull = 25
    # Constraints based on problem statement:
    # v1, v2, v3, v_bull are distinct positive integers.
    # v1 < v2 < v3 < v_bull
    # s1 <= 35  => v1 <= 7
    # s_bull >= 50 => v_bull >= 10
    # Bobby's total / 5 = 46. Cliff's total / 5 = 37. Number of shots = 5.
    
    possible_bull_scores = set()
    
    # Iterate over possible v1 values (v1 from 1 to 7)
    for v1 in range(1, 8):
        # Iterate over possible v2 values.
        # From Anna's eq: v_bull = 25 - 3*v1 - v2
        # Since v_bull >= 10 => 25 - 3*v1 - v2 >= 10 => v2 <= 15 - 3*v1
        v2_max = 15 - 3 * v1
        for v2 in range(v1 + 1, v2_max + 1):
            v_bull = 25 - 3 * v1 - v2
            
            # The existence of v3 requires v_bull >= v2 + 2
            if v_bull < v2 + 2:
                continue
                
            # Iterate through possible v3 values to find if there exists one
            # that makes the score set valid.
            for v3 in range(v2 + 1, v_bull):
                v_scores = (v1, v2, v3, v_bull)
                
                # Check if Cliff's score (37) can be formed.
                global memo
                memo.clear()
                cliff_possible = can_form_sum(37, v_scores, 5)
                
                if not cliff_possible:
                    continue
                
                # Check if Bobby's score (46) can be formed.
                memo.clear()
                bobby_possible = can_form_sum(46, v_scores, 5)
                
                if bobby_possible:
                    # A valid score set is found. Add the bull's eye score to the set.
                    # s_bull = v_bull * 5
                    possible_bull_scores.add(v_bull * 5)
                    # Break from the v3 loop as we only need one working v3
                    # to prove this s_bull is possible.
                    break
    
    sorted_scores = sorted(list(possible_bull_scores))
    equation_parts = ["1"] * len(sorted_scores)
    
    print("Based on the provided information, the possible values for the bull's eye score are:")
    print(sorted_scores)
    print("\nThe final equation is the sum of each possibility to find the total count:")
    print(f"{' + '.join(equation_parts)} = {len(sorted_scores)}")
    print(f"\nThe exact number of possible values for the score of the bull's eye is {len(sorted_scores)}.")
    
solve()
<<<11>>>