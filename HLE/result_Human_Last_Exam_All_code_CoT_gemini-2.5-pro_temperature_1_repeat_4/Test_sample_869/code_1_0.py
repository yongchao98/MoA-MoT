import math
from itertools import permutations

# A global list to store all valid distribution matrices
matrices = []

def factorial(n):
    """Computes factorial."""
    if n == 0:
        return 1
    return math.factorial(n)

def generate_matrices_recursive(matrix, row, col, row_sums, col_sums):
    """
    Recursively generates all 5x5 matrices with given row and column sums.
    """
    n = 5
    if row == n:
        # If all rows are filled, check if the matrix is valid
        matrices.append([r[:] for r in matrix])
        return

    # Calculate the next position
    next_row, next_col = (row, col + 1) if col + 1 < n else (row + 1, 0)

    # Determine the range of possible values for the current cell
    # The value cannot exceed the remaining sum for the current row or column
    max_val = min(row_sums[row], col_sums[col])
    
    # Pruning for the last column
    if col == n - 1:
        if row_sums[row] == col_sums[col]:
             val = row_sums[row]
             matrix[row][col] = val
             row_sums[row] -= val
             col_sums[col] -= val
             generate_matrices_recursive(matrix, next_row, next_col, row_sums, col_sums)
             row_sums[row] += val
             col_sums[col] += val
    else:
        for val in range(max_val + 1):
            matrix[row][col] = val
            row_sums[row] -= val
            col_sums[col] -= val
            generate_matrices_recursive(matrix, next_row, next_col, row_sums, col_sums)
            # Backtrack
            row_sums[row] += val
            col_sums[col] += val

def is_favorable(c):
    """
    Checks if a distribution matrix C is favorable.
    A matrix is favorable if there's a permutation p such that for each
    person i, they are the champion for type p[i].
    """
    n = 5
    # Iterate through all possible champion-type assignments (permutations)
    for p in permutations(range(n)):
        is_p_favorable = True
        # For the current permutation, check if every person is a champion
        for i in range(n):
            champion_type = p[i]
            person_i_count = c[i][champion_type]
            
            is_champion = True
            for k in range(n):
                if i == k:
                    continue
                if person_i_count <= c[k][champion_type]:
                    is_champion = False
                    break
            
            if not is_champion:
                is_p_favorable = False
                break
        
        if is_p_favorable:
            return True # Found a valid permutation, so the matrix is favorable
            
    return False

def count_sequences(c):
    """
    Calculates the number of sequences N(C) for a given matrix C.
    N(C) = (5!^5) / product(c_ij!)
    """
    n = 5
    prod_cij_fact = 1
    for i in range(n):
        for j in range(n):
            prod_cij_fact *= factorial(c[i][j])
    
    # Using integer division
    return factorial(n)**n // prod_cij_fact

def solve_probability():
    """
    Main function to solve the problem.
    """
    n = 5
    target_sum = 5
    
    # 1. Generate all valid distribution matrices
    initial_matrix = [[0] * n for _ in range(n)]
    row_sums = [target_sum] * n
    col_sums = [target_sum] * n
    generate_matrices_recursive(initial_matrix, 0, 0, row_sums, col_sums)

    # 2. Calculate S, the total number of sequences
    s_numerator = factorial(n * target_sum)
    s_denominator = factorial(target_sum)**n
    S = s_numerator // s_denominator

    # 3. Calculate F, the number of favorable sequences
    F = 0
    for c in matrices:
        if is_favorable(c):
            F += count_sequences(c)

    # 4. Simplify and print the result
    common_divisor = math.gcd(F, S)
    f_simple = F // common_divisor
    s_simple = S // common_divisor
    
    print(f"Total number of possible distributions (S) is calculated as 25! / (5!)^5.")
    print(f"S = {s_numerator} / {s_denominator} = {S}")
    print(f"The number of favorable distributions (F) is the sum of N(C) for all favorable matrices C.")
    print(f"F = {F}")
    print(f"The probability P = F / S.")
    print(f"P = {F} / {S} = {f_simple} / {s_simple}")

solve_probability()
