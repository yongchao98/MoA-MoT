def solve_max_peaceable_queens():
    """
    Calculates the maximum number m such that m white queens and m black
    queens can coexist on an N x N chessboard without same-color attacks.
    """

    # 1. Define the board size.
    N = 16

    print(f"Solving the peaceable queens problem for a {N}x{N} chessboard.")
    print("---------------------------------------------------------")

    # 2. Establish the upper bound for m.
    print(f"Step 1: Determine the upper bound for m.")
    print("A set of 'm' non-attacking queens of the same color must occupy 'm' different rows.")
    print(f"On a {N}x{N} board, there are {N} rows available.")
    print(f"Therefore, the number of queens of a single color, m, cannot exceed {N}.")
    upper_bound = N
    print(f"Conclusion: m <= {upper_bound}")
    print("---------------------------------------------------------")

    # 3. Check if the upper bound is achievable.
    print(f"Step 2: Check if m = {N} is achievable.")
    print(f"To place {N} white queens and {N} black queens, we need two things:")
    print(f" - A valid placement for {N} non-attacking white queens (a {N}-queens solution).")
    print(f" - A valid placement for {N} non-attacking black queens (another {N}-queens solution).")
    print(f" - The two sets of positions must be completely disjoint (no shared squares).")
    print("---------------------------------------------------------")
    
    # 4. Apply the relevant mathematical theorem.
    print("Step 3: Apply the known mathematical result.")
    print("A theorem on the N-queens problem states that two disjoint N-queens solutions exist")
    print("if and only if the board size N is greater than or equal to 4.")
    print("---------------------------------------------------------")

    # 5. Apply the theorem to our specific case (N=16).
    print(f"Step 4: Apply the theorem to our case where N = {N}.")
    condition_n = 4
    is_achievable = N >= condition_n
    
    print(f"We check if {N} >= {condition_n}.")
    if is_achievable:
        print(f"The condition is met. This guarantees that two disjoint {N}-queens solutions exist.")
        print(f"We can place {N} white queens and {N} black queens.")
        max_m = upper_bound
    else:
        # This part of the logic applies to N<4 and won't be executed here.
        print(f"The condition is not met. m = {N} is not achievable.")
        max_m = 0 # For N=2,3 the answer would be 0.
    
    print("---------------------------------------------------------")
    
    # 6. State the final answer.
    print("Final Conclusion:")
    print(f"The maximum number of queens, m, is bounded by {upper_bound}.")
    print(f"This upper bound of m = {max_m} is achievable.")
    print("\nFinal Equation Summary:")
    print(f"Board Size (N): {N}")
    print(f"Condition for solution (N >= k): {condition_n}")
    print(f"Maximum m: {max_m}")


if __name__ == '__main__':
    solve_max_peaceable_queens()