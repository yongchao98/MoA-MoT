import math

def calculate_minimum_prototypes(N):
    """
    Calculates the minimum total number of prototypes for N+1 concentric circles.

    The problem requires perfectly classifying N+1 concentric circles (classes 0 to N)
    with a 1-NN classifier using a balanced number of prototypes per class.

    Args:
        N (int): The index of the largest circle. The number of circles is N+1.
                 N must be a non-negative integer.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    # For N=0, there is only one class (a point at the origin). 
    # One prototype is sufficient. The formula handles this case correctly.
    # k = floor(pi*0)+1 = 1. Total = (0+1)*1 = 1.
    
    # Step 1: Calculate the minimum number of prototypes per class, k.
    # This is derived from the geometric condition: k > pi * N
    k = math.floor(math.pi * N) + 1
    
    # Step 2: Calculate the total number of classes.
    # The circles are indexed t=0, 1, ..., N.
    num_classes = N + 1
    
    # Step 3: Calculate the total number of prototypes.
    # This is (number of classes) * (prototypes per class).
    total_prototypes = num_classes * k
    
    # Print the breakdown of the calculation
    print(f"For N = {N}:")
    print(f"Number of classes = {N} + 1 = {num_classes}")
    print(f"Prototypes per class (k) = floor(pi * {N}) + 1 = {k}")
    print(f"Total number of prototypes = ({N} + 1) * {k} = {total_prototypes}")
    
    # Return the final answer in the required format
    print(f"\n<<<{total_prototypes}>>>")


# --- User Input ---
# You can change the value of N here to solve for different numbers of circles.
# For example, if you have 3 circles (radii 0, c, 2c), then N=2.
N_value = 2
calculate_minimum_prototypes(N_value)
