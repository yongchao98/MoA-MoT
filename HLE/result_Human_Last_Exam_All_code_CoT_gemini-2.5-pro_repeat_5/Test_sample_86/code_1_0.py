import math

def calculate_min_prototypes(N):
    """
    Calculates the minimum total number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.

    Args:
        N (int): The number of circles excluding the center point (t=1, 2, ..., N).
                 The total number of classes is N+1.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # For N=0, there is 1 class (the origin). We need 1 prototype to classify it.
    # The formula also works for N=0: k = floor(0) + 1 = 1. Total = (0+1)*1 = 1.
    
    # Step 1: Calculate the minimum number of prototypes per class (k).
    # This is derived from the geometric constraint to ensure perfect classification.
    # k must satisfy: k > pi * N
    # The smallest integer k is floor(pi * N) + 1.
    k = math.floor(math.pi * N) + 1

    # Step 2: Calculate the total number of classes.
    # The circles are indexed t=0, 1, ..., N, so there are N+1 classes.
    num_classes = N + 1

    # Step 3: Calculate the total minimum number of prototypes.
    total_prototypes = num_classes * k

    # Step 4: Output the results, showing the numbers in the final equation.
    print(f"For a dataset with N = {N}:")
    print(f"The total number of classes is N + 1 = {num_classes}.")
    print(f"The minimum required prototypes per class (k) is floor(π * N) + 1 = {k}.")
    print("-" * 30)
    print(f"The minimum total number of prototypes is:")
    print(f"({num_classes}) * ({k}) = {total_prototypes}")


# --- Example Usage ---
# You can change the value of N here to solve for a different number of circles.
N_value = 12
calculate_min_prototypes(N_value)