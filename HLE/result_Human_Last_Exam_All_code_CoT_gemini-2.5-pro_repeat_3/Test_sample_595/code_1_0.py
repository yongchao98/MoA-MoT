import math

def solve():
    """
    Finds the maximum number of grid squares a triangle with sides 18, 18, 18*sqrt(2)
    can pass through without its perimeter containing any lattice points.
    """
    
    # Case 1: Legs are aligned with the grid axes.
    # This corresponds to vectors (18, 0) and (0, 18).
    # Leg 1 (length 18 on an axis) crosses 18 grid lines.
    # Leg 2 (length 18 on an axis) crosses 18 grid lines.
    # The hypotenuse from (18,0) to (0,18) crosses 18 vertical and 18 horizontal lines.
    # Total k is the sum of all crossings.
    k_axis_aligned = 18 + 18 + 18 + 18
    max_k = k_axis_aligned
    best_A, best_B = 18, 0

    # Case 2: General orientation.
    # The number of squares k = 4*floor(a) + 2*floor(b) + 1, where a^2+b^2 = 18^2.
    # We search for integers A = floor(a) and B = floor(b) that maximize k.
    # We can assume a > b, so A >= B.
    # The maximum value for a is less than 18, so A can be at most 17.
    for A in range(18):
        # From A = floor(a) and a^2+b^2=324, we have A <= a < A+1.
        # This gives a constraint on B = floor(b).
        # A^2 + B^2 <= a^2 + b^2 = 324.
        if 324 - A*A < 0:
            continue
        B_max_val = math.sqrt(324 - A*A)
        
        # We only need to check B <= A.
        for B in range(A + 1):
            if B > B_max_val:
                continue

            # We need to ensure that the circle a^2+b^2=324 actually intersects
            # the square [A, A+1) x [B, B+1).
            # A sufficient condition is that the circle does not lie completely
            # outside or completely inside the larger square [A, B] to [A+1, B+1].
            # Closest point in square to origin is (A,B). A^2+B^2 <= 324 is already checked.
            # Furthest point is (A+1,B+1). We need (A+1)^2+(B+1)^2 > 324.
            if (A + 1)**2 + (B + 1)**2 <= 324:
                continue

            # Calculate k for this valid (A, B) pair.
            current_k = 4 * A + 2 * B + 1
            if current_k > max_k:
                max_k = current_k
                best_A = A
                best_B = B

    print("The number of squares crossed, k, is maximized when the triangle is placed in a general orientation.")
    print("The formula for k is 4*A + 2*B + 1, where A=floor(a) and B=floor(b) for the leg vector (a,b).")
    print(f"The search finds the optimal integer parts to be A = {best_A} and B = {best_B}.")
    print(f"The final equation for the maximum k is:")
    print(f"k = 4 * {best_A} + 2 * {best_B} + 1 = {max_k}")

solve()
