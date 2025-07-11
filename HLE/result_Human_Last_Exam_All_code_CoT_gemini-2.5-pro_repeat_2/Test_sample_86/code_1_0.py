import math

def calculate_prototypes():
    """
    Calculates the minimum number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.
    """
    try:
        # Prompt the user to enter the value of N
        N_str = input("Enter the value of N (number of circles is N+1): ")
        N = int(N_str)
        if N < 0:
            print("Error: N must be a non-negative integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer for N.")
        return

    # There are N+1 classes in total.
    num_classes = N + 1

    # Handle the edge case where N=0.
    # We have one class (a single point at the origin).
    # We need at least one prototype per class, so k=1.
    if N == 0:
        k = 1
        total_prototypes = 1
        print("\nFor N = 0, we have 1 class (the origin).")
        print("The minimum number of prototypes per class (k) required is: 1")
        print(f"The total number of prototypes is {k} * {num_classes} = {total_prototypes}")
        return

    # For N >= 1, we must satisfy the condition derived from the geometry
    # of the two outermost circles (N-1 and N).
    # The condition is: cos(pi/k) > (2N - 1) / (2N)
    
    # Calculate the value on the right side of the inequality
    val = (2 * N - 1) / (2 * N)

    # To solve for k, we rearrange the inequality:
    # pi/k < arccos((2N - 1) / (2N))
    # k > pi / arccos((2N - 1) / (2N))
    
    # Calculate the angle whose cosine is val
    angle = math.acos(val)
    
    # Calculate the minimum k as a float
    k_float = math.pi / angle
    
    # Since k must be an integer and the inequality is strict,
    # k must be the next integer greater than k_float.
    # We use floor(k_float) + 1 to find this integer.
    # This correctly handles cases where k_float is an integer,
    # as k must be strictly greater.
    k = math.floor(k_float) + 1
    
    # The total number of prototypes is k times the number of classes.
    total_prototypes = k * num_classes

    # Print the results, showing each number in the final equation
    print(f"\nFor N = {N}, we have {num_classes} classes.")
    print(f"The minimum number of prototypes per class (k) required is: {k}")
    print(f"The total number of prototypes is {k} * {num_classes} = {total_prototypes}")

if __name__ == '__main__':
    calculate_prototypes()