import math
import argparse

def solve_prototypes(N):
    """
    Calculates the minimum total number of prototypes for N+1 concentric circles.
    
    Args:
        N (int): The parameter defining the number of circles (from t=0 to N).
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # Handle the base case N=0 separately for clarity.
    if N == 0:
        print("For N=0, there is only one class (a single point at the center).")
        print("Therefore, only 1 prototype is needed.")
        print("\nFinal Answer for N=0:")
        print("Total = (0 + 1) * 1 = 1")
        return

    # Let k be the number of prototypes per class. The derived condition for
    # perfect classification is 2*N*sin(pi/(2k)) < 1.
    # This leads to the minimum integer k being floor(pi * N) + 1.
    
    # Calculate k, the number of prototypes per class
    pi_times_N = math.pi * N
    k = math.floor(pi_times_N) + 1
    
    # Calculate the total number of prototypes
    num_classes = N + 1
    total_prototypes = num_classes * k

    print(f"To guarantee perfect classification for N = {N}:")
    print("-" * 40)
    
    print("Step 1: Calculate the minimum prototypes required per class (k).")
    print("The formula is: k = floor(π * N) + 1")
    print(f"k = floor({math.pi:.6f} * {N}) + 1")
    print(f"k = floor({pi_times_N:.6f}) + 1")
    print(f"k = {math.floor(pi_times_N)} + 1 = {k}\n")
    
    print("Step 2: Calculate the total number of prototypes.")
    print("The formula is: Total = (Number of Classes) * k = (N + 1) * k")
    # As requested, showing each number in the final equation
    print(f"Total = ({N} + 1) * {k} = {num_classes} * {k} = {total_prototypes}")


def main():
    """Main function to parse arguments and run the solver."""
    parser = argparse.ArgumentParser(
        description="Calculate the minimum number of prototypes for classifying N+1 concentric circles with a 1-NN classifier.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'N', 
        type=int, 
        help="A non-negative integer specifying the outermost circle's index 't'.\nThere will be N+1 circles in total, from t=0 to N."
    )
    args = parser.parse_args()
    solve_prototypes(args.N)


if __name__ == "__main__":
    main()