import cmath

def solve_and_sum():
    """
    This function solves the given system of matrix equations and finds the sum
    of the first coordinates of all solutions.
    """
    # The first equation is:
    # | 5  0 | X1^2 + X1^2 | 6  0 | = | -53/12  0 |
    # | 0 -5 |           | 0  6 |   |     0   0 |
    #
    # This simplifies to (A + 6I) * X1^2 = C1.
    # Let the first coordinate of X1 be x1. Then the (1,1) element equation is:
    # (5 + 6) * x1^2 = -53/12
    # 11 * x1^2 = -53/12
    # x1^2 = -53 / (12 * 11)
    x1_squared = -53 / 132

    # The second equation is:
    # | 4  0 | X2^2 + X2^2 | 6  0 | = | -3/11  0 |
    # | 0 -5 |           | 0  6 |   |    0   0 |
    #
    # This simplifies to (A' + 6I) * X2^2 = C2.
    # Let the first coordinate of X2 be x2. Then the (1,1) element equation is:
    # (4 + 6) * x2^2 = -3/11
    # 10 * x2^2 = -3/11
    # x2^2 = -3 / (11 * 10)
    x2_squared = -3 / 110

    # The solutions for the first coordinate of X1 are the square roots of x1_squared.
    x1_sol1 = cmath.sqrt(x1_squared)
    x1_sol2 = -x1_sol1

    # The solutions for the first coordinate of X2 are the square roots of x2_squared.
    x2_sol1 = cmath.sqrt(x2_squared)
    x2_sol2 = -x2_sol1

    # The set of all possible first coordinates are {x1_sol1, x1_sol2, x2_sol1, x2_sol2}.
    all_solutions = [x1_sol1, x1_sol2, x2_sol1, x2_sol2]
    total_sum = sum(all_solutions)

    # Print the sum equation with each number.
    print("The first equation gives solutions for the first coordinate of X1.")
    print(f"The first solution for X1's first coordinate is: {x1_sol1}")
    print(f"The second solution for X1's first coordinate is: {x1_sol2}")
    print("\nThe second equation gives solutions for the first coordinate of X2.")
    print(f"The first solution for X2's first coordinate is: {x2_sol1}")
    print(f"The second solution for X2's first coordinate is: {x2_sol2}")
    
    print("\nThe sum of all these first coordinates is:")
    print(f"({x1_sol1}) + ({x1_sol2}) + ({x2_sol1}) + ({x2_sol2}) = {total_sum}")

    # The sum is expected to be 0 since we are summing pairs of opposite numbers.
    # We take the real part to present a clean final answer.
    final_answer = total_sum.real
    print(f"\nFinal numerical sum: {final_answer}")

if __name__ == "__main__":
    solve_and_sum()