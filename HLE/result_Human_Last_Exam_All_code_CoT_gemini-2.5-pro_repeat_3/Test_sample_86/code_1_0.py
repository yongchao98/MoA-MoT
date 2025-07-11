import math

def calculate_min_prototypes(N):
    """
    Calculates the minimum number of prototypes for N+1 concentric circles.

    Args:
        N (int): The number of circles is N+1, with radii tc for t=0,1,...,N.
                 N must be a non-negative integer.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # Case N=0: Only one class (the origin).
    # The minimum number of balanced prototypes is k=1 per class.
    if N == 0:
        k = 1
        num_classes = 1
        total_prototypes = 1
        print(f"For N = {N}, there is only {num_classes} class (a point at the origin).")
        print(f"The minimum number of prototypes per class (k) is {k}.")
        print(f"Final Equation: {k} * ({N} + 1) = {total_prototypes}")
        return

    # Case N > 0:
    # The number of prototypes per class, k, must satisfy k > pi / acos(1 - 1/(2*N)).
    # We calculate this lower bound for k.
    try:
        min_k_float = math.pi / math.acos(1.0 - 1.0 / (2.0 * N))
    except ValueError:
        # This can happen if the argument to acos is outside [-1, 1],
        # though it shouldn't for N > 0.
        print("Error during calculation. Please check the value of N.")
        return

    # k must be the smallest integer strictly greater than this value.
    k = math.floor(min_k_float) + 1

    # The total number of prototypes is k multiplied by the number of classes (N+1).
    num_classes = N + 1
    total_prototypes = k * num_classes

    print(f"For a dataset with N = {N}, there are {num_classes} concentric circles (classes).")
    print("To guarantee perfect classification with a 1-NN classifier and balanced classes,")
    print(f"the number of prototypes per class, k, must be greater than {min_k_float:.4f}.")
    print(f"The smallest integer k that satisfies this is {k}.")
    print("The minimum total number of prototypes is k * (N + 1).")
    print(f"Final Equation: {k} * {num_classes} = {total_prototypes}")


# --- Main execution ---
# You can change the value of N here to test different scenarios.
N_value = 10
calculate_min_prototypes(N_value)