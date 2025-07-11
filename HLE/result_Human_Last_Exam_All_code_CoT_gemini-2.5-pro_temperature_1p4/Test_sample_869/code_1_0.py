import math
from itertools import permutations
import sys

def main():
    """
    This script calculates the probability described in the problem.
    The problem involves calculating two large numbers:
    S: The total number of ways to distribute the items.
    F: The number of 'favorable' ways, where each person has a unique 'dominant' item type.
    The final probability is P = F / S.
    
    The script will print the final equation showing the values for F and S.
    """
    
    # Increase recursion limit for the deep search required to find F.
    if sys.getrecursionlimit() < 2000:
        sys.setrecursionlimit(2000)

    # --- Step 1: Calculate S ---
    # S represents the total number of ways to distribute 25 items (5 items each of 5 types)
    # among 5 individuals, where each individual receives 5 items. This is calculated as 25! / (5!)^5.
    
    S = math.factorial(25) // (math.factorial(5)**5)

    # --- Step 2: Calculate F ---
    # F is the number of favorable distributions. A distribution is favorable if for each person,
    # there is a unique type of item they hold more of than any other person.
    # We model a distribution by a 5x5 matrix C, where C[i][j] is the count of type j items for person i.
    
    # We find all valid matrices C for a fixed dominance assignment (Person i dominates Type i),
    # calculate the number of shuffles (Ways(C)) for each, and sum them up. This sum is F_fixed.
    # The total F is then 5! * F_fixed due to symmetry.
    # Ways(C) = (5!)^5 / product_over_all_i_j(C[i][j]!)

    # Pre-calculate partitions and permutations for rows to build matrices
    partitions = [
        (5,0,0,0,0), (4,1,0,0,0), (3,2,0,0,0), (3,1,1,0,0), 
        (2,2,1,0,0), (2,1,1,1,0), (1,1,1,1,1)
    ]
    
    all_rows = []
    for p in partitions:
        all_rows.extend(list(set(permutations(p))))

    fact_cache = {i: math.factorial(i) for i in range(6)}

    # This global variable will accumulate the ways for all valid matrices.
    F_fixed_ways = 0
    
    memo_ways = {}
    def calculate_ways(matrix_tuple):
        """Calculates the number of shuffles for a given distribution matrix."""
        if matrix_tuple in memo_ways:
            return memo_ways[matrix_tuple]
        
        denominator = 1
        for r in range(5):
            for c in range(5):
                denominator *= fact_cache[matrix_tuple[r][c]]
        
        res = (fact_cache[5]**5) // denominator
        memo_ways[matrix_tuple] = res
        return res

    def check_dominance(matrix):
        """Checks if the diagonal dominance condition holds for a matrix."""
        for j in range(5):  # For each column (type)
            diag_val = matrix[j][j]
            for i in range(5):  # For each row (individual)
                if i == j: continue
                if matrix[i][j] >= diag_val:
                    return False
        return True

    def find_favorable_matrices(k, matrix, col_sums):
        """Recursively builds and checks matrices."""
        nonlocal F_fixed_ways
        
        # Base case: matrix is fully constructed
        if k == 5:
            if check_dominance(matrix):
                matrix_tuple = tuple(map(tuple, matrix))
                F_fixed_ways += calculate_ways(matrix_tuple)
            return

        # Optimization: The last row is determined by the previous rows' column sums.
        if k == 4:
            last_row = tuple(5 - cs for cs in col_sums)
            if any(v < 0 for v in last_row) or sum(last_row) != 5:
                return 
            
            find_favorable_matrices(k + 1, matrix + [list(last_row)], (5,5,5,5,5))
            return

        # Recursive step for rows 0 through 3
        for row in all_rows:
            new_col_sums = tuple(cs + r for cs, r in zip(col_sums, row))
            if all(cs <= 5 for cs in new_col_sums):
                find_favorable_matrices(k + 1, matrix + [list(row)], new_col_sums)

    # --- Execute the Search ---
    print("Calculating F... (Note: This is computationally intensive and may take a few minutes)")
    find_favorable_matrices(0, [], (0,0,0,0,0))
    print("Calculation of F complete.")
    
    # Final calculation of F
    F = math.factorial(5) * F_fixed_ways

    # --- Final Output ---
    # We will now print each number in the final equation P = F/S
    print("\n--- Results ---")
    print(f"The total number of ways to distribute the items is:")
    print(f"S = {S}")
    
    print(f"\nThe number of favorable distributions is:")
    print(f"F = {F}")
    
    print(f"\nThe probability P = F / S is:")
    
    common_divisor = math.gcd(F, S)
    F_simplified = F // common_divisor
    S_simplified = S // common_divisor
    
    print(f"P = {F} / {S} = {F_simplified} / {S_simplified}")
    
    final_prob_val = F/S
    print(f"P ≈ {final_prob_val:.10f}")
    
if __name__ == '__main__':
    main()