import sys

# It's recommended to run this with n and m as command-line arguments.
# Example: python your_script_name.py 3 3
# If not provided, it will use default values n=3, m=3.

def count_rank_k_matrices(rows, cols, k):
    """
    Calculates the number of `rows x cols` matrices of rank `k` over GF(2).
    """
    if k < 0 or k > min(rows, cols):
        return 0
    if k == 0:
        return 1

    # Numerator part for rows
    num_rows = 1
    for i in range(k):
        num_rows *= (2**rows - 2**i)

    # Numerator part for cols
    num_cols = 1
    for i in range(k):
        num_cols *= (2**cols - 2**i)

    # Denominator part
    den = 1
    for i in range(k):
        den *= (2**k - 2**i)

    return (num_rows * num_cols) // den

def count_p_positions(n, m):
    """
    Calculates the number of losing positions (P-positions) for an n x m grid.
    """
    if n == 0 or m == 0:
        return 1
    
    total_p_positions = 0
    rows = n - 1
    cols = m - 1
    
    limit = min(rows, cols)
    
    for k in range(limit + 1):
        num_rank_k = count_rank_k_matrices(rows, cols, k)
        total_p_positions += num_rank_k * (4**k)
        
    return total_p_positions

def main():
    """
    Main function to calculate and print the result.
    """
    if len(sys.argv) == 3:
        try:
            n = int(sys.argv[1])
            m = int(sys.argv[2])
        except ValueError:
            print("Invalid input. Please provide two integers for n and m.")
            return
    else:
        # Default values if no command-line arguments are given
        n, m = 3, 3

    if n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return

    num_p_positions = count_p_positions(n, m)
    
    # The condition for the first player to have a winning strategy with >50% probability
    # is that the number of P-positions is less than half the total number of matrices.
    threshold = 2**(n * m - 1)
    
    # The problem asks to output the numbers in the final equation.
    # The final equation is the inequality comparing num_p_positions and threshold.
    is_winning_prob_gt_50 = num_p_positions < threshold
    
    print(f"For an {n}x{m} matrix:")
    print(f"The number of P-positions is |S_P| = {num_p_positions}")
    print(f"The threshold is 2^(n*m - 1) = {threshold}")
    
    comparison_char = '<' if is_winning_prob_gt_50 else '>='
    print(f"The condition for the first player to have a high winning probability is |S_P| < 2^(n*m-1).")
    print(f"Checking the inequality: {num_p_positions} {comparison_char} {threshold}")

    if is_winning_prob_gt_50:
        print("The probability is strictly more than 50%. The function f(n, m) returns 1.")
    else:
        print("The probability is not strictly more than 50%. The function f(n, m) returns 0.")


if __name__ == "__main__":
    main()
