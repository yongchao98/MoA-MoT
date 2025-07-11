def is_valid_solution(queens, n):
    """
    Checks if a given set of queen positions is a valid N-Queens solution.
    A solution is valid if no two queens attack each other.
    The coordinates are 1-based (1 to N).
    """
    if len(queens) != n:
        return False
    
    rows = set()
    cols = set()
    diag1 = set()  # r - c
    diag2 = set()  # r + c

    for r, c in queens:
        if not (1 <= r <= n and 1 <= c <= n):
             return False # Check if queens are on the board
        if r in rows or c in cols or (r - c) in diag1 or (r + c) in diag2:
            return False  # Found an attack
        rows.add(r)
        cols.add(c)
        diag1.add(r - c)
        diag2.add(r + c)
        
    return True

def construct_solution_set(n):
    """
    Constructs a solution for the N-Queens problem for specific N.
    This construction works for N where N % 6 is not 2 or 3.
    Since 16 % 6 = 4, this construction is valid.
    Coordinates are 1-based.
    """
    queens = set()
    # Place queens in the first half of rows
    for i in range(1, n // 2 + 1):
        queens.add((i, 2 * i))
    # Place queens in the second half of rows
    for i in range(1, n // 2 + 1):
        queens.add((n // 2 + i, 2 * i - 1))
    return queens

def reflect_solution_vertically(queens, n):
    """
    Reflects a given solution across the vertical midline of the board.
    A queen at (r, c) moves to (r, n + 1 - c).
    """
    reflected_queens = set()
    for r, c in queens:
        reflected_queens.add((r, n + 1 - c))
    return reflected_queens

def solve_problem():
    """
    Main function to solve the problem and print the reasoning.
    """
    n = 16
    
    print(f"Problem: Find the maximum number m such that m white and m black queens can coexist on a {n}x{n} board.")
    print("This is interpreted as queens of the same color not attacking each other.\n")

    # Step 1: Determine the theoretical maximum for m
    print("Step 1: Determine the upper bound for m.")
    print(f"The maximum number of non-attacking queens of one color on a {n}x{n} board is {n}.")
    print(f"Therefore, m cannot be greater than {n}. The maximum possible value for m is {n}.\n")

    # Step 2: Prove that m = 16 is achievable
    print("Step 2: Prove that m = 16 is achievable by construction.")
    print("We will construct two disjoint sets of 16 non-attacking queen positions.\n")

    # Construction
    white_queens = construct_solution_set(n)
    black_queens = reflect_solution_vertically(white_queens, n)
    print("Constructed a set of positions for 16 white queens (S1).")
    print("Constructed a second set for 16 black queens (S2) by reflecting S1 vertically.\n")

    # Step 3: Verification
    print("Step 3: Verify the constructions.")
    valid_white = is_valid_solution(white_queens, n)
    print(f"  - Verification 1: Are the {len(white_queens)} white queens in S1 non-attacking? Result: {valid_white}")
    
    valid_black = is_valid_solution(black_queens, n)
    print(f"  - Verification 2: Are the {len(black_queens)} black queens in S2 non-attacking? Result: {valid_black}")
    
    are_disjoint = white_queens.isdisjoint(black_queens)
    print(f"  - Verification 3: Are S1 and S2 disjoint (no shared squares)? Result: {are_disjoint}\n")

    # Step 4: Conclusion
    print("Step 4: Final Conclusion.")
    if valid_white and valid_black and are_disjoint:
        m = n
        print(f"We have successfully found two disjoint non-attacking sets of {n} queens.")
        print(f"This proves that m = {n} is achievable.")
        print(f"Since the maximum possible value for m is {n}, the answer must be {n}.\n")
        print("Final Answer Equation:")
        print(f"max_m = {m}")
    else:
        print("The construction and verification failed. The answer could not be determined with this method.")

if __name__ == "__main__":
    solve_problem()