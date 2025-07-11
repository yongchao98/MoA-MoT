import math

def solve():
    """
    Calculates the minimum total number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.
    """
    try:
        n_str = input("Enter the value of N (the label of the outermost circle, N>=0): ")
        N = int(n_str)
        if N < 0:
            print("N must be a non-negative integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Handle the special case for N=0
    if N == 0:
        # If N=0, there is only one class (t=0, the origin).
        # One prototype is sufficient for perfect classification.
        k = 1
        num_classes = 1
        total_prototypes = 1
    else:
        # For N >= 1, we must satisfy the condition for all adjacent classes.
        # The most restrictive condition comes from the boundary between
        # class N-1 and class N.
        # The condition is: cos(pi/k) > (N - 0.5) / N
        
        # We need to find the smallest integer k such that the inequality holds.
        # k > pi / acos((N - 0.5) / N)
        
        # Right-hand side of the inequality
        rhs = (N - 0.5) / N
        
        # Solve for k
        # We need pi/k < acos(rhs), so k > pi / acos(rhs)
        k_float = math.pi / math.acos(rhs)
        
        # k must be an integer, so we take the floor and add 1.
        k = math.floor(k_float) + 1
        
        # Total number of prototypes is k * (number of classes)
        num_classes = N + 1
        total_prototypes = k * num_classes

    print("\n--- Calculation ---")
    print(f"For N = {N}, there are {N+1} classes (from 0 to {N}).")
    print(f"The minimum number of balanced prototypes required per class (k) is: {k}")
    print(f"The number of classes (N+1) is: {N+1}")
    print(f"The final equation for the total number of prototypes is: M = k * (N+1)")
    print(f"So, M = {k} * {N+1} = {total_prototypes}")
    
    # Final answer in the required format
    print(f"\nMinimum total number of prototypes:")
    print(f'<<<{total_prototypes}>>>')

solve()