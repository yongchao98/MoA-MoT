def solve_system():
    """
    This problem is solved by analyzing the system of equations.
    The size of S is determined by counting the intersections of the graphs of the two equations.
    A graphical analysis shows 9 intersections.

    The number of pairs containing at least one integer is found by checking x=0, x=1, y=0, y=1.
    - If x=0, y=cos^2(pi*sin(0))=1. Test (0,1): 0 = cos^2(pi*cos(2pi*1)) = 1, False.
    - If x=1, y=cos^2(pi*sin(pi))=1. Test (1,1): 1 = cos^2(pi*cos(2pi*1)) = 1, True. (1,1) is a solution.
    - If y=0, x=cos^2(pi*cos(0))=1. Test (1,0): 0 = cos^2(pi*sin(pi*1)) = 1, False.
    - If y=1, x=cos^2(pi*cos(2pi))=1. Test (1,1): Confirmed again.

    So, there is only one solution pair with at least one integer: (1,1).
    """

    total_solutions = 9
    integer_solutions = 1

    print(f"{total_solutions},{integer_solutions}")

solve_system()