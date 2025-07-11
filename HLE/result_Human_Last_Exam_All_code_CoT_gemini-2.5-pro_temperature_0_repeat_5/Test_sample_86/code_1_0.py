import math

def solve_prototypes():
    """
    Calculates the minimum number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.
    """
    # The user can change this value.
    # N is the number of circles with radius > 0.
    # The total number of classes is N + 1.
    N = 5

    print(f"Solving for N = {N}")
    print("-" * 20)

    # Handle the edge case where N=0
    if N == 0:
        # If N=0, there is only one class (the origin).
        # We need at least one prototype. The "balanced" constraint is trivial.
        num_classes = 1
        k = 1
        total_prototypes = 1
        print(f"For N = {N}, there is {num_classes} class.")
        print(f"The minimum number of prototypes required is 1.")
        return

    # For N > 0
    num_classes = N + 1

    # Calculate the minimum number of prototypes per class, k.
    # The condition is: k > pi / arccos(1 - 1/(2*N))
    try:
        arccos_val = math.acos(1 - 1 / (2 * N))
        # k must be an integer, so we take the floor and add 1.
        k = math.floor(math.pi / arccos_val) + 1
    except ValueError:
        print("Error: Invalid value for arccos. N must be a positive number.")
        return
    except ZeroDivisionError:
        print("Error: N cannot be zero in this formula part (handled by the N=0 case).")
        return


    # Calculate the total number of prototypes
    total_prototypes = num_classes * k

    print(f"Number of classes (N+1): {num_classes}")
    print(f"Minimum prototypes per class (k): {k}")
    print("\nThe final equation for the total number of prototypes is:")
    print(f"Total = (N+1) * k")
    print(f"Total = {num_classes} * {k}")
    print(f"Result: {total_prototypes}")

solve_prototypes()