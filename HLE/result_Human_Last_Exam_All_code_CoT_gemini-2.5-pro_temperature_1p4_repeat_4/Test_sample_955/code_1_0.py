import math

def solve_for_n():
    """
    This function finds the value of n based on the problem description.
    """
    
    # The probability equation is derived as:
    # (Number of reachable odd cells + Number of reachable even cells) / Total cells = 0.66
    # (n^2/2 + 5n - 12) / n^2 = 0.66
    # This simplifies to the quadratic equation: 0.16n^2 - 5n + 12 = 0
    
    # Solve the quadratic equation 0.16n^2 - 5n + 12 = 0
    a, b, c = 0.16, -5, 12
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        # This case should not be reached based on the values
        print("The equation has no real solutions.")
        return

    # Calculate the two possible roots for n
    n1 = (-b + math.sqrt(discriminant)) / (2*a)
    n2 = (-b - math.sqrt(discriminant)) / (2*a)

    # The problem implies n is a reasonably large, even integer. We take the larger root.
    n_solution = n1 if n1 > n2 else n2

    # Since the root is not an integer, we test the nearest even integers.
    n_floor = math.floor(n_solution)
    n_ceil = math.ceil(n_solution)

    # Find the two even integers bracketing the solution
    n_low = n_floor if n_floor % 2 == 0 else n_floor - 1
    n_high = n_ceil if n_ceil % 2 == 0 else n_ceil + 1
    
    # Define a function to calculate the probability for a given n
    def probability(n):
        if n == 0:
            return 0
        reachable_odds = (n**2) / 2
        # This formula is valid for even n >= 4
        reachable_evens = 5*n - 12
        total_reachable = reachable_odds + reachable_evens
        return total_reachable / (n**2)

    prob_low = probability(n_low)
    prob_high = probability(n_high)
    
    target_prob = 0.66
    
    # Determine which of the two even integers gives a probability closer to 66%
    if abs(prob_low - target_prob) < abs(prob_high - target_prob):
        final_n = n_low
    else:
        final_n = n_high
        
    # Output the final result
    num_odd = int((final_n**2) / 2)
    num_even = 5 * final_n - 12
    total_reachable_cells = num_odd + num_even
    total_grid_cells = final_n**2
    final_prob = total_reachable_cells / total_grid_cells
    
    print("The final equation is:")
    print(f"(({final_n}^2 / 2) + (5 * {final_n} - 12)) / {final_n}^2 = {final_prob:.4f}")
    print(f"({num_odd} + {num_even}) / {total_grid_cells} = {final_prob:.4f}")
    print(f"{total_reachable_cells} / {total_grid_cells} = {final_prob:.4f}")
    print("\nThis probability is the closest to the given 66%.")
    print(f"The determined value of n is: {final_n}")


solve_for_n()