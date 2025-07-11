import math

def solve_for_n():
    """
    Solves for the grid size 'n' based on the given probability.
    
    The relationship is R / (n*n) = P, where R is the number of reachable cells,
    n is the grid size, and P is the probability.
    
    We are given:
    - P = 0.66
    - n must be an even integer.
    - R must be an integer (as it's a count of cells).

    This means we are looking for the smallest even integer 'n' for which
    0.66 * n^2 is a whole number.
    """
    
    probability = 0.66
    found_n = 0
    
    # We start searching from n=2, checking only even numbers.
    n_candidate = 2
    while True:
        # Calculate the corresponding number of reachable cells, R.
        R = probability * (n_candidate ** 2)
        
        # Check if R is an integer. We use a small tolerance for floating-point inaccuracies.
        if abs(R - round(R)) < 1e-9:
            found_n = n_candidate
            break
        
        # Move to the next even number.
        n_candidate += 2

    # Calculate the final integer value for R.
    final_R = int(round(probability * (found_n ** 2)))
    
    # As requested, output each number in the final equation: R / (n * n) = P
    print(f"The calculated number of reachable cells (R) is: {final_R}")
    print(f"The determined grid size (n) is: {found_n}")
    print("\nBased on these values, the equation is:")
    print(f"{final_R} / ({found_n} * {found_n}) = {probability}")

solve_for_n()
<<<10>>>