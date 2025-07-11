import math

def solve_grid_problem():
    """
    Solves for the grid size 'n' based on the given probability.

    The problem can be modeled with the following steps:
    1. The total number of cells in the grid is n*n.
    2. The number of "white" cells (where x+y is odd) is n*n / 2. All of these are reachable within 2 moves.
    3. The number of reachable "black" cells (where x+y is even) is calculated to be 5*n - 12 for even n >= 6.
    4. The total probability P(n) is the sum of reachable white and black cells, divided by the total cells:
       P(n) = ( (n*n / 2) + (5*n - 12) ) / (n*n)
    5. We set this probability to the given value of 66% (0.66) and solve for n.
       (n*n/2 + 5*n - 12) / (n*n) = 0.66
       0.5 + (5*n - 12) / (n*n) = 0.66
       (5*n - 12) / (n*n) = 0.16
       5*n - 12 = 0.16 * n*n
       This gives the quadratic equation: 0.16*n^2 - 5*n + 12 = 0
    """

    # Coefficients for the quadratic equation: a*n^2 + b*n + c = 0
    a = 0.16
    b = -5
    c = 12

    # Calculate the discriminant
    discriminant = b**2 - 4*a*c

    # Since the problem implies a solution exists, but the quadratic does not yield an integer,
    # we find the closest even integer 'n' that satisfies the condition.
    best_n = -1
    min_diff = float('inf')
    target_prob = 0.66

    # We test a range of even integers for n.
    # The quadratic solution is ~28.63, so we check around this value.
    for n in range(6, 50, 2):
        reachable_white = (n * n) / 2
        # This formula is valid for even n >= 6
        reachable_black = 5 * n - 12
        total_reachable = reachable_white + reachable_black
        
        current_prob = total_reachable / (n * n)
        diff = abs(current_prob - target_prob)

        if diff < min_diff:
            min_diff = diff
            best_n = n

    # Final calculation with the best found 'n'
    n = best_n
    reachable_white_final = int((n * n) / 2)
    reachable_black_final = 5 * n - 12
    total_reachable_final = reachable_white_final + reachable_black_final
    total_cells_final = n * n
    
    print("The problem is to find n for the equation:")
    print(f"( (n*n / 2) + (5*n - 12) ) / (n*n) = {target_prob}")
    print("\nSince solving this for n does not yield an integer, we find the even integer n that gives the closest probability.")
    print(f"Testing n = {n}:")
    print(f"Reachable White Cells = {n}*{n}/2 = {reachable_white_final}")
    print(f"Reachable Black Cells = 5*{n} - 12 = {reachable_black_final}")
    print(f"Total Reachable Cells = {reachable_white_final} + {reachable_black_final} = {total_reachable_final}")
    print(f"Total Grid Cells = {n}*{n} = {total_cells_final}")
    print(f"\nThe full probability equation for n={n} is:")
    print(f"(({reachable_white_final}) + ({reachable_black_final})) / ({total_cells_final}) = {total_reachable_final/total_cells_final:.4f}")
    print(f"\nThis probability is the closest to {target_prob}.")
    print(f"The value of n is {n}.")


solve_grid_problem()
<<<28>>>