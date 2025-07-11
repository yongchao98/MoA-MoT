import sys
from itertools import combinations_with_replacement
from collections import Counter

def solve_archery_puzzle():
    """
    Solves the archery puzzle to find the number of possible values for the bull's eye score.
    """

    # Bobby's and Cliff's total scores, divided by 5 for simplified calculation
    k_bobby_total = 230 // 5
    k_cliff_total = 185 // 5

    # There are 56 ways to distribute 5 arrows into 4 bins (target areas).
    # We pre-calculate these partitions to avoid re-computing them.
    # Each partition is a tuple (n1, n2, n3, n_bull) where the sum is 5.
    partitions = []
    for combo in combinations_with_replacement(range(4), 5):
        counts = Counter(combo)
        partitions.append(tuple(counts.get(i, 0) for i in range(4)))

    def is_achievable(k_list, k_total):
        """
        Checks if a k_total can be achieved with 5 arrows given a list of score multipliers.
        k_list: A tuple (k1, k2, k3, k_bull)
        k_total: The target sum of multipliers (e.g., 46 for Bobby).
        """
        k1, k2, k3, k_bull = k_list
        for p in partitions:
            n1, n2, n3, n_bull = p
            current_sum = n1 * k1 + n2 * k2 + n3 * k3 + n_bull * k_bull
            if current_sum == k_total:
                return True
        return False

    valid_bull_eye_scenarios = {}

    # From Anna's equation (3*k1 + k2 + k_bull = 25) and constraints (0 < k1 < k2 < k_bull),
    # we can deduce that k_bull must be in the range [8, 20].
    for k_bull in range(8, 21):
        if k_bull in valid_bull_eye_scenarios:
            continue
        
        # From Anna's equation: k2 = 25 - k_bull - 3*k1
        # Iterate through possible k1 values.
        for k1 in range(1, k_bull):
            k2 = 25 - k_bull - (3 * k1)
            
            # Check if this k2 satisfies the strict inequality k1 < k2.
            if k1 < k2:
                # Also, there must be space for at least one integer k3 between k2 and k_bull.
                # So, k_bull must be at least k2 + 2.
                if k2 < k_bull - 1:
                    # We have a valid (k1, k2, k_bull) triplet from Anna's data.
                    # Now, we must find if there exists a k3 that makes Bobby's and Cliff's scores possible.
                    for k3 in range(k2 + 1, k_bull):
                        k_list = (k1, k2, k3, k_bull)
                        
                        bobby_ok = is_achievable(k_list, k_bobby_total)
                        if bobby_ok:
                            cliff_ok = is_achievable(k_list, k_cliff_total)
                            if cliff_ok:
                                # Found a valid scenario for this k_bull.
                                # Store it and break the inner loops to check the next k_bull.
                                s_bull = 5 * k_bull
                                s1 = 5 * k1
                                s2 = 5 * k2
                                valid_bull_eye_scenarios[s_bull] = (s1, s2)
                                break  # out of k3 loop
            if (5 * k_bull) in valid_bull_eye_scenarios:
                break # out of k1 loop

    # Print the final results
    num_solutions = len(valid_bull_eye_scenarios)
    print(f"The exact number of possible values for the score of the bull's eye is {num_solutions}.")
    
    if num_solutions > 0:
        print("\nFor each possible bull's eye score, one valid example for Anna's score is:")
        sorted_s_bulls = sorted(valid_bull_eye_scenarios.keys())
        for s_bull in sorted_s_bulls:
            s1, s2 = valid_bull_eye_scenarios[s_bull]
            total = 3 * s1 + s2 + s_bull
            print(f"Bull's Eye = {s_bull}: 3 * {s1} + 1 * {s2} + 1 * {s_bull} = {total}")

    sys.stdout.flush() # Ensure all output is printed
    return num_solutions


if __name__ == '__main__':
    final_answer = solve_archery_puzzle()
    print(f"\n<<<{final_answer}>>>")
