import math

def calculate_min_prototypes(N):
    """
    Calculates the minimum total number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.

    Args:
        N (int): The number of circles is N+1, with radii tc for t=0,1,...,N.

    Returns:
        None. Prints the result.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # Handle the special case for N=0, where there is only one class (a point).
    # One prototype at the origin is sufficient.
    if N == 0:
        num_classes = 1
        prototypes_per_class = 1
        total_prototypes = 1
        print(f"For N=0, we have {num_classes} class.")
        print(f"The number of prototypes per class is {prototypes_per_class}.")
        print(f"The minimum total number of prototypes is: {num_classes} * {prototypes_per_class} = {total_prototypes}")
        return

    # For N > 0, we derive the minimum number of prototypes per class, k.
    # The condition for perfect separation is: sin(pi / (2*k)) < 1 / (2 * sqrt(N))
    # This means k must be > pi / (2 * asin(1 / (2 * sqrt(N))))
    
    # Calculate the value for the arcsin function
    val = 1 / (2 * math.sqrt(N))
    
    # Calculate the critical value for k (where equality holds)
    k_critical = math.pi / (2 * math.asin(val))
    
    # The number of prototypes per class, k, must be the smallest integer
    # strictly greater than k_critical.
    prototypes_per_class = math.floor(k_critical) + 1
    
    # The number of classes is N + 1 (for t = 0, 1, ..., N)
    num_classes = N + 1
    
    # Total number of prototypes is k * (N + 1)
    total_prototypes = prototypes_per_class * num_classes
    
    print(f"For N={N}, we have {num_classes} classes.")
    print(f"The minimum number of prototypes required per class is {prototypes_per_class}.")
    print("The minimum total number of prototypes is:")
    print(f"{num_classes} (classes) * {prototypes_per_class} (prototypes per class) = {total_prototypes}")


# --- User Input ---
# You can change the value of N here to see the result for different numbers of circles.
N_value = 10
# --- End of User Input ---

calculate_min_prototypes(N_value)