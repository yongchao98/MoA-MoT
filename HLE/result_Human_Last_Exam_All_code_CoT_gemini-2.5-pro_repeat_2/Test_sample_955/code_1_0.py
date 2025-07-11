import math

def solve_for_n():
    """
    This function solves for the grid size 'n' based on the problem description.
    """
    target_prob = 0.66

    print("Step 1: Deriving the formula for the number of reachable cells, N(n).")
    print("The total number of reachable cells N(n) for an n x n grid (where n is an even integer >= 6) is given by the formula N(n) = 6n - 5.")
    print("This is derived by counting all unique cells reachable in at most 3 moves:")
    print(" - White cells reachable: (n+2) from initial diagonals + 5 from D-B-B paths = n + 7")
    print(" - Black cells reachable: 5n - 12 (from diagonals passing through the 7 D-B path endpoints)")
    print(" - Total N(n) = (n + 7) + (5n - 12) = 6n - 5.\n")
    
    print(f"Step 2: Set up the probability equation N(n) / n^2 = {target_prob}.")
    print(f"(6n - 5) / n^2 = {target_prob}")
    print(f"This gives the quadratic equation: {target_prob}n^2 - 6n + 5 = 0\n")

    print("Step 3: Solve the quadratic equation.")
    a, b, c = 0.66, -6, 5
    discriminant = b**2 - 4*a*c
    if discriminant >= 0:
        n1 = (-b + math.sqrt(discriminant)) / (2*a)
        n2 = (-b - math.sqrt(discriminant)) / (2*a)
        print(f"The solutions are n = {n1:.2f} and n = {n2:.2f}.")
    print("Since the solutions are not integers, the 66% probability is likely an approximation.\n")

    print("Step 4: Test even integer values of n to find the closest probability.")
    
    best_n = 0
    min_diff = float('inf')

    # We test for even n starting from 4.
    # The formula 6n-5 is valid for n>=6. n=4 is a special case.
    
    # Case n=4
    n = 4
    # For n=4, all 16 cells are reachable.
    reachable_cells_4 = 16
    prob_4 = reachable_cells_4 / (n**2)
    diff = abs(prob_4 - target_prob)
    print(f"For n = {n}: Reachable cells = {reachable_cells_4}, Total = {n**2}, P({n}) = {prob_4:.4f}, |P(n) - 0.66| = {diff:.4f}")
    if diff < min_diff:
        min_diff = diff
        best_n = n

    # Cases n >= 6
    for n in range(6, 21, 2):
        reachable_cells = 6 * n - 5
        total_cells = n**2
        prob = reachable_cells / total_cells
        diff = abs(prob - target_prob)
        print(f"For n = {n}: Reachable cells = {reachable_cells}, Total = {total_cells}, P({n}) = {prob:.4f}, |P(n) - 0.66| = {diff:.4f}")
        if diff < min_diff:
            min_diff = diff
            best_n = n

    print(f"\nThe value of n that results in a probability closest to {target_prob} is n = {best_n}.\n")
    
    print("Step 5: Final probability equation for the determined value of n.")
    final_n = best_n
    final_reachable = 6 * final_n - 5
    final_total = final_n**2
    final_prob = final_reachable / final_total
    
    print(f"For n = {final_n}, the final equation is:")
    print(f"{final_reachable} / {final_total} = {final_prob:.6f}")
    
    return final_n

result_n = solve_for_n()
print(f"<<<{result_n}>>>")