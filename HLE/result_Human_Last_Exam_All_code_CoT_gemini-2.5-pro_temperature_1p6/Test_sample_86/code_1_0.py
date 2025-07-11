import math

def calculate_min_prototypes(N):
    """
    Calculates the minimum total number of balanced prototypes for a 1-NN classifier
    to perfectly classify N+1 concentric circles.

    Args:
        N (int): The number of circles, where circle radii are c, 2c, ..., Nc.
                 The total number of classes is N+1 (including the center point).

    Returns:
        tuple: A tuple containing k (prototypes per class) and the total number of prototypes.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return None, None

    # The central point (t=0) counts as one of the N+1 classes.
    # If N=0, there is only one class (the origin). One prototype is sufficient.
    # The balanced class rule requires k prototypes. k must be at least 1.
    # Our formula works, as pi*sqrt(0)=0 -> floor(0)+1=1, max(4,1)=4. 
    # But for N=0, 1 prototype is enough. k=1 is a valid interpretation for k>pi*sqrt(0).
    # However, k>3 is a general constraint for any circle to be separated from the center.
    # If N=0 there are no circles to separate from the center. 1 is sufficient.
    if N == 0:
        k = 1
        total_prototypes = 1
        return k, total_prototypes

    # Constraint from separating any circle from the origin: k > 3
    # Constraint from separating adjacent circles: k > pi * sqrt(N)
    # We take the maximum of these two constraints.
    # k must be an integer.
    k = math.floor(math.pi * math.sqrt(N)) + 1
    
    # The actual k must also be at least 4 due to the origin separation constraint
    # (unless N=0).
    k = max(4, k)
    
    num_classes = N + 1
    total_prototypes = num_classes * k
    
    return k, total_prototypes

# You can change the value of N here
N = 10

# Perform the calculation
k_per_class, total = calculate_min_prototypes(N)

# Print the result showing each number in the equation
if k_per_class is not None:
    print(f"For N = {N}:")
    print(f"The number of classes is N + 1 = {N + 1}.")
    print(f"The minimum number of prototypes required per class (k) is {k_per_class}.")
    print(f"The minimum total number of prototypes is (N + 1) * k = {N + 1} * {k_per_class} = {total}.")
    # Final answer format for N=10
    if N == 10:
        print("\n<<<110>>>")