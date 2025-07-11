import math

def find_grid_size():
    """
    This function determines the grid size 'n' based on the problem's probability condition.

    It first calculates the non-integer solution to the derived probability equation
    to find candidate values for n. Then, it tests the nearest even integers to find
    the one that satisfies the condition that the probability is 66% (rounded).
    Finally, it prints the step-by-step verification for the correct 'n'.
    """

    # The derived probability equation is 0.16*n^2 - 5*n + 12 = 0.
    # We solve this to find the approximate value of n.
    a = 0.16
    b = -5
    c = 12
    
    # Calculate the discriminant
    delta = b**2 - 4 * a * c
    
    # The problem implies a real solution must exist
    if delta < 0:
        print("The equation has no real solutions.")
        return

    # Calculate the potential values for n using the quadratic formula
    n_solution = (-b + math.sqrt(delta)) / (2 * a)

    # Since n must be an even integer, we test the even integers closest to the solution.
    # The relevant solution is the larger one (~28.63).
    # The closest even integers are 28 and 30.
    candidate_values = [28, 30]
    
    best_n = None
    
    for n_test in candidate_values:
        # The derived formula for the number of reachable cells (R) is n^2/2 + 5n - 12.
        # This formula is valid for even n >= 6.
        if n_test < 6 or n_test % 2 != 0:
            continue
        
        reachable_cells = (n_test**2) / 2 + 5 * n_test - 12
        total_cells = n_test**2
        probability = reachable_cells / total_cells
        
        # Check if the calculated probability rounds to 0.66
        if round(probability, 2) == 0.66:
            best_n = n_test
            break # Found the solution
            
    if best_n is None:
        print("Could not find an even integer n that satisfies the condition.")
        return

    # Final calculation and output for the determined value of n
    n = best_n
    reachable_cells = int((n**2) / 2 + 5 * n - 12)
    total_cells = int(n**2)
    final_probability = reachable_cells / total_cells

    print(f"The grid size is determined to be n x n, where n = {n}.")
    print("This was found by deriving the formula for the number of reachable cells (R) and solving the equation R/n^2 ≈ 0.66.")
    print("-" * 30)
    print("Final Verification:")
    print(f"For n = {n}:")
    print(f"1. The number of reachable 'odd' cells is n^2 / 2 = {n**2 // 2}.")
    print(f"2. The number of reachable 'even' cells is 5*n - 12 = 5*{n} - 12 = {5*n - 12}.")
    print(f"3. The total number of reachable cells (R) is {n**2 // 2} + {5*n - 12} = {reachable_cells}.")
    print(f"4. The total number of cells in the grid (N) is n^2 = {n}^2 = {total_cells}.")
    print("\nThe final equation using these numbers is:")
    print(f"Probability = Reachable Cells / Total Cells")
    print(f"Probability = {reachable_cells} / {total_cells} = {final_probability:.4f}")
    print(f"\nThis probability, {final_probability:.4f}, rounds to 0.66, which matches the problem statement.")
    
# Execute the function
find_grid_size()
<<<28>>>