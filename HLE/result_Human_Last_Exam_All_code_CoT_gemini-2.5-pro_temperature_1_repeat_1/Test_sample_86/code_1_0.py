import math

def solve_prototypes():
    """
    Calculates the minimum total number of prototypes for a given N.
    """
    # Set the number of outer circles, N.
    # For the problem given, N+1 is the total number of circles (classes).
    # Please change the value of N here to solve for different cases.
    N = 5

    if N < 0:
        print("N must be a non-negative integer.")
        return

    print(f"For N = {N}:")

    # Handle the special case where N=0 (a single point-class)
    if N == 0:
        num_classes = 1
        k = 1
        total_prototypes = 1
        print(f"Number of classes = {num_classes}")
        print(f"Minimum prototypes per class (k) = {k}")
        print(f"Total minimum prototypes = {k} * {num_classes} = {total_prototypes}")
        return

    # For N >= 1
    num_classes = N + 1

    # The condition for k is: cos(pi/k) > 1 - 1/(2*N)
    # This is equivalent to k > pi / acos(1 - 1/(2*N))
    # We need to find the smallest integer k that satisfies this.

    # The value inside the arccos
    val = 1.0 - (1.0 / (2.0 * N))
    
    # Calculate the threshold for k
    arccos_val = math.acos(val)
    k_float = math.pi / arccos_val
    
    # k must be an integer strictly greater than k_float
    k = math.floor(k_float) + 1

    total_prototypes = k * num_classes

    print(f"Number of classes = N + 1 = {num_classes}")
    print(f"Minimum prototypes per class (k) = {k}")
    print(f"Total minimum prototypes = k * (N + 1) = {k} * {num_classes} = {total_prototypes}")

# Execute the function
solve_prototypes()