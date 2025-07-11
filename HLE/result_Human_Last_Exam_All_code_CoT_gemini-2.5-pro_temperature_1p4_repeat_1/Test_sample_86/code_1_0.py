import math

def solve_prototypes():
    """
    Calculates the minimum number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.
    """
    # Set the value for N, the index of the outermost circle.
    # The total number of circles (and classes) is N + 1.
    N = 5

    print(f"Solving for N = {N}")
    print("-------------------")

    if N < 0:
        print("N must be a non-negative integer.")
        return

    # Handle the trivial case where there is only one class (the origin).
    if N == 0:
        k = 1
        total_prototypes = 1
        print("For N=0, there is only one class (the origin).")
        print("Minimum prototypes per class (k) is 1.")
        print(f"Total number of prototypes = {k} * ({N} + 1) = {total_prototypes}")
        return

    # For N > 0, we use the derived formula to find k, the number of prototypes per class.
    # The condition is: cos(pi / k) > 1 - 1 / (2N)
    # This is equivalent to: k > pi / acos(1 - 1 / (2N))

    # 1. Calculate the threshold value
    threshold = 1.0 - 1.0 / (2.0 * N)
    print(f"The condition for the number of prototypes per class (k) is: cos(pi / k) > 1 - 1/(2*{N}) = {threshold:.4f}")

    # 2. Solve for the lower bound of k
    # This value represents the exact (non-integer) k where the decision boundary would touch the points.
    k_bound = math.pi / math.acos(threshold)
    print(f"Solving for k gives the inequality: k > pi / acos({threshold:.4f}) ≈ {k_bound:.4f}")

    # 3. Find the minimum integer k
    # k must be an integer strictly greater than k_bound.
    k = math.floor(k_bound) + 1
    print(f"The smallest integer k satisfying this inequality is {k}.")

    # 4. Calculate the total number of prototypes
    # Total = k * (Number of classes) = k * (N + 1)
    num_classes = N + 1
    total_prototypes = k * num_classes

    print("\nFinal Calculation:")
    print(f"Minimum prototypes per class (k): {k}")
    print(f"Total number of classes (N+1): {num_classes}")
    print(f"Minimum total number of prototypes = {k} * {num_classes} = {total_prototypes}")


# Run the solver
solve_prototypes()