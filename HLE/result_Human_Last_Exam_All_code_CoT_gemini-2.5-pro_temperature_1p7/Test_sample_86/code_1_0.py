import math

# The user can change this value to calculate the result for a different N.
N = 20

def calculate_min_prototypes(N):
    """
    Calculates the minimum total number of prototypes needed to perfectly
    classify N+1 concentric circles with a 1-NN classifier using a balanced
    set of prototypes lying on the circles.
    
    Args:
        N: The number of circles excluding the center point (t=0). 
           N must be a non-negative integer.
    """
    print(f"Solving for N = {N}")
    print("-" * 30)

    if not isinstance(N, int) or N < 0:
        print("Error: N must be a non-negative integer.")
        return None

    # Handle the special case where N=0.
    # We have one class (a single point at the origin), so one prototype is sufficient.
    if N == 0:
        k = 1
        total_prototypes = 1
        print("For N=0, we have one class (a point at the origin).")
        print("We need k=1 prototype per class.")
        print("\nFinal Equation:")
        print(f"Total Prototypes = (N + 1) * k")
        print(f"Total Prototypes = ({N} + 1) * {k} = {total_prototypes}")
        return total_prototypes

    # For N >= 1, the number of prototypes per class, k, is determined by the
    # condition that for any point on any circle t, its nearest prototype is of class t.
    # The most restrictive condition arises from ensuring that points on the outermost
    # circle (t=N) are not misclassified as belonging to the adjacent inner circle (t=N-1).
    # This leads to the inequality: cos(pi / k) > 1 - 1 / (2*N)
    
    print("The minimum number of prototypes per class (k) is found from the inequality:")
    print("cos(pi / k) > 1 - 1 / (2 * N)")
    
    # We can solve this inequality for k:
    # pi / k < arccos(1 - 1 / (2 * N))
    # k > pi / arccos(1 - 1 / (2 * N))
    
    print("\nStep-by-step calculation to find the smallest integer k:")
    
    # 1. Calculate the right-hand side of the inequality
    rhs = 1.0 - 1.0 / (2.0 * N)
    print(f"  1. cos(pi / k) > 1 - 1 / (2 * {N})  =>  cos(pi / k) > {rhs:.6f}")
    
    # 2. Take arccos of the value.
    #    This gives the upper bound for the angle pi/k.
    arc_cos_val = math.acos(rhs)
    print(f"  2. pi / k < arccos({rhs:.6f})  =>  pi / k < {arc_cos_val:.6f} radians")
    
    # 3. Solve for the lower bound of k
    k_threshold = math.pi / arc_cos_val
    print(f"  3. k > pi / {arc_cos_val:.6f}  =>  k > {k_threshold:.6f}")
    
    # 4. The minimum k is the smallest integer greater than this threshold.
    k = math.floor(k_threshold) + 1
    print(f"\n==> The smallest integer k must be {k}.")
    
    # Calculate the total number of prototypes required
    num_classes = N + 1
    total_prototypes = num_classes * k
    
    print("\nFinal calculation for the total number of prototypes:")
    print(f"Number of classes = N + 1 = {num_classes}")
    print(f"Prototypes per class = k = {k}")
    print("\nFinal Equation:")
    print(f"Total Prototypes = (N + 1) * k")
    print(f"Total Prototypes = ({N} + 1) * {k} = {total_prototypes}")
    
    return total_prototypes

# --- Main execution ---
final_answer = calculate_min_prototypes(N)
if final_answer is not None:
    print(f"\nFor N = {N}, the minimum total number of prototypes required is {final_answer}.")
    print(f"<<<{final_answer}>>>")
