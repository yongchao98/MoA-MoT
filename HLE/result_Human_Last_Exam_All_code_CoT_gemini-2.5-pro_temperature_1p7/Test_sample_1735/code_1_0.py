import itertools
import math

def calculate_max_cumulative(seq):
    """Calculates the maximum absolute cumulative sum for a sequence."""
    if not seq:
        return 0, [], []
    current_sum = 0
    cumulative_sums = []
    abs_cumulative_sums = []
    for x in seq:
        current_sum += x
        cumulative_sums.append(current_sum)
        abs_cumulative_sums.append(abs(current_sum))
    return max(abs_cumulative_sums) if abs_cumulative_sums else 0, cumulative_sums, abs_cumulative_sums

def solve_greedy_f(A):
    """
    Implements the greedy strategy from statement (F):
    At each step k, choose x_k to minimize |C(k)|.
    """
    remaining = list(A)
    greedy_seq = []
    current_sum = 0
    
    while remaining:
        best_next_num = None
        min_abs_sum = float('inf')
        
        for num in remaining:
            abs_sum = abs(current_sum + num)
            if abs_sum < min_abs_sum:
                min_abs_sum = abs_sum
                best_next_num = num
        
        greedy_seq.append(best_next_num)
        current_sum += best_next_num
        remaining.remove(best_next_num)
        
    return greedy_seq

def solve_optimal_bruteforce(A):
    """Finds the optimal permutation by checking all unique permutations."""
    min_max_c = float('inf')
    optimal_seq = []
    
    # Use set to handle duplicate permutations if input has duplicate numbers
    for p in set(itertools.permutations(A)):
        max_c, _, _ = calculate_max_cumulative(p)
        if max_c < min_max_c:
            min_max_c = max_c
            optimal_seq = list(p)
            
    return optimal_seq

def main():
    """
    Demonstrates that the greedy algorithm from statement (F) is not optimal.
    """
    # A = [3, -1, -4, 2]
    # A = [1, 1, -2]
    # A = [10, -10, 1, -1] # Counterexample for greedy approach F
    A = [1, -1, 10, -10]

    print(f"Analyzing sequence A = {A}\n")
    
    # --- Part 1: Run the greedy algorithm from statement (F) ---
    greedy_seq = solve_greedy_f(A)
    max_c_greedy, s_greedy, c_greedy = calculate_max_cumulative(greedy_seq)
    
    print("--- Verifying Statement (F) ---")
    print("Greedy approach: Minimize |C(k)| at each step")
    print(f"Greedy Permutation: {greedy_seq}")
    print(f"Cumulative Sums: {s_greedy}")
    print(f"Absolute Cumulative Sums: {c_greedy}")
    print(f"Resulting max{{C(k)}}: {max_c_greedy}\n")
    
    # --- Part 2: Find the true optimal solution via brute-force ---
    optimal_seq = solve_optimal_bruteforce(A)
    max_c_optimal, s_optimal, c_optimal = calculate_max_cumulative(optimal_seq)

    print("--- Brute-Force Optimal Solution ---")
    print(f"Optimal Permutation: {optimal_seq}")
    print(f"Cumulative Sums: {s_optimal}")
    print(f"Absolute Cumulative Sums: {c_optimal}")
    print(f"Optimal max{{C(k)}}: {max_c_optimal}\n")
    
    # --- Part 3: Conclusion ---
    if max_c_greedy > max_c_optimal:
        print("Conclusion: The greedy approach from (F) is not optimal, as shown above.")
    else:
        print("Conclusion: For this specific case, the greedy approach was optimal.")

if __name__ == "__main__":
    main()
