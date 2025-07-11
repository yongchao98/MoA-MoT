import sympy

def solve_bounds():
    """
    This function determines the lower and upper bounds for the real number t
    based on the problem's mathematical requirements.

    The logic is as follows:
    1.  Let A be the interval of possible sums a_i + a_j, which is [-2, 2*t].
    2.  The core condition is that for any x in A, 1/x must also be in A.
    3.  This first requires that 0 is not in A, which implies t < 0.
    4.  The condition then becomes that the image of A under f(x)=1/x is a subset of A.
    5.  The image of A = [-2, 2*t] is f(A) = [1/(2*t), -1/2].
    6.  The requirement f(A) subset A leads to the inequalities:
        - 1/(2*t) >= -2  => t <= -1/4
        - -1/2 <= 2*t    => t >= -1/4
    7.  These inequalities force t to be exactly -1/4.
    """

    # The only value for t that satisfies the conditions is -1/4.
    t_solution = -1/4

    # The lower and upper bounds of a set containing a single point are the point itself.
    lower_bound = t_solution
    upper_bound = t_solution

    # The problem asks to output numbers from the "final equation".
    # The final equation from the problem statement is (a_0 + a_2)(a_1 + a_3) = 1.
    # The numeric constants defining the problem are 1 (from the equation)
    # and -1 (from the interval [-1, t]).
    # We print these as requested.
    print("The numeric constants in the problem statement are 1 and -1.")
    print(f"The derived value for t is {t_solution}.")

    # Print the lower and upper bounds separated by a space.
    print(f"{lower_bound} {upper_bound}")

solve_bounds()