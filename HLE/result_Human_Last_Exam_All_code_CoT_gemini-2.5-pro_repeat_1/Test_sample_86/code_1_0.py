import math

def calculate_min_prototypes(N):
    """
    Calculates the minimum total number of prototypes to perfectly classify
    N+1 concentric circles with a 1-NN classifier under a balanced class constraint.

    Args:
        N (int): The number of circles is N+1, with labels 0, 1, ..., N.
                 N must be a non-negative integer.
    """
    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return

    print(f"Solving for N = {N}")
    print("-" * 20)

    # Handle the base case: N=0 means only one class (a point at the origin).
    # One prototype is sufficient.
    if N == 0:
        num_classes = 1
        k = 1
        total_prototypes = 1
        print("For N=0, there is only one class (a single point at the origin).")
        print("A single prototype is sufficient.")
        print(f"The final equation is: {num_classes} * {k} = {total_prototypes}")
        print(f"\nMinimum total number of prototypes needed: {total_prototypes}")
        return

    # For N > 0, we derive the condition for the number of prototypes per class, k.
    # The condition for perfect separation of any two adjacent circles t-1 and t is:
    # cos(pi / k) > 1 - 1 / (2*t)
    # This condition is most stringent for the largest t, i.e., t = N.
    num_classes = N + 1
    print(f"Number of classes = N + 1 = {num_classes}")
    print("The condition for the number of prototypes per class (k) is:")
    print("cos(pi / k) > 1 - 1 / (2 * N)")

    # Calculate the right-hand side of the inequality
    rhs = 1 - 1 / (2 * N)
    print(f"For N = {N}, the condition is: cos(pi / k) > 1 - 1 / {2 * N}")
    print(f"cos(pi / k) > {rhs}")

    # To solve for k, we rearrange the inequality:
    # pi / k < arccos(rhs)
    # k > pi / arccos(rhs)
    # The arccos is well-defined since 0.5 <= rhs < 1 for N >= 1.
    angle = math.acos(rhs)
    k_threshold = math.pi / angle
    
    print("\nSolving for k:")
    print(f"k > pi / arccos({rhs:.4f})")
    print(f"k > {math.pi:.4f} / {angle:.4f}")
    print(f"k > {k_threshold:.4f}")

    # k must be the smallest integer greater than this threshold.
    k = math.floor(k_threshold) + 1
    print(f"The minimum integer k satisfying this is {k}.")

    # The total number of prototypes is (N+1) * k
    total_prototypes = num_classes * k

    print("\nCalculating the total number of prototypes:")
    print(f"Total prototypes = (N + 1) * k")
    # Output the final equation with the calculated numbers
    print("The final equation is:")
    print(f"({N} + 1) * {k} = {total_prototypes}")
    
    print(f"\nMinimum total number of prototypes needed: {total_prototypes}")


# Example usage: Set the value of N here.
# The user can change this value to solve for a different N.
N_val = 10
calculate_min_prototypes(N_val)

# The general formula for the answer is (N+1) * (floor(pi / arccos(1 - 1/(2N))) + 1) for N>0, and 1 for N=0.
# For N=10, the answer is 110.
# For N=1, k = floor(pi/acos(0.5)) + 1 = floor(3)+1 = 4. Total = (1+1)*4=8.
# For N=2, k = floor(pi/acos(0.75)) + 1 = floor(4.34)+1 = 5. Total = (2+1)*5=15.
# For N=20, k = floor(pi/acos(0.975)) + 1 = floor(13.99)+1 = 14. Total = (20+1)*14=294.

# Let's use the derived result to provide the final answer
N = 10
k = math.floor(math.pi / math.acos(1 - 1/(2*N))) + 1
total = (N+1) * k
print(f'<<<{total}>>>')