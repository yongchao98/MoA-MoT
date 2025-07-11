from fractions import Fraction

def solve_variance():
    """
    Calculates the variance of Y, where Y is the second closest value to X1
    among X2, X3, and X4, interpreted as the second smallest distance.
    """
    # Step 1 & 2: Define the problem and the distribution of distances.
    # Let D_i = |X_{i+1} - X_1|. Y is the second order statistic of {D_1, D_2, D_3}.
    # For D = |U1 - U2| with U1, U2 ~ U[0,1], the CDF is F_D(d) = 2d - d^2
    # and the PDF is f_D(d) = 2 - 2d for d in [0, 1].

    # Step 3: The PDF of Y = D_(2) is f_Y(y) = 6 * F_D(y) * (1 - F_D(y)) * f_D(y).
    # f_Y(y) = 6 * (2y - y^2) * (1 - (2y - y^2)) * (2 - 2y)
    # f_Y(y) = 12y * (2 - y) * (1 - y)^3

    # Step 4: Calculate E[Y] and E[Y^2] by integrating y*f_Y(y) and y^2*f_Y(y).
    # These integrals can be solved analytically, for example using the Beta function.
    # E[Y] = integral from 0 to 1 of 12y^2(2-y)(1-y)^3 dy = 11/35
    # E[Y^2] = integral from 0 to 1 of 12y^3(2-y)(1-y)^3 dy = 9/70

    E_Y_num = 11
    E_Y_den = 35
    E_Y = Fraction(E_Y_num, E_Y_den)

    E_Y_squared_num = 9
    E_Y_squared_den = 70
    E_Y_squared = Fraction(E_Y_squared_num, E_Y_squared_den)

    # Step 5: Calculate the variance using Var(Y) = E[Y^2] - (E[Y])^2.
    variance = E_Y_squared - E_Y**2

    # Step 6: Print the results, showing each number in the final equation.
    print("This script calculates the variance of Y, interpreted as the second smallest distance.")
    print("The formula for the variance is Var(Y) = E[Y^2] - (E[Y])^2.")
    print(f"The expected value, E[Y], is calculated to be {E_Y_num}/{E_Y_den}.")
    print(f"The expected value of Y squared, E[Y^2], is calculated to be {E_Y_squared_num}/{E_Y_squared_den}.")
    
    print("\nCalculating the variance:")
    # The prompt requires printing each number in the final equation.
    print(f"Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"Var(Y) = {E_Y_squared_num}/{E_Y_squared_den} - ({E_Y_num}/{E_Y_den})^2")
    print(f"Var(Y) = {E_Y_squared.numerator}/{E_Y_squared.denominator} - {E_Y.numerator**2}/{E_Y.denominator**2}")
    print(f"Var(Y) = {variance.numerator}/{variance.denominator}")

    final_answer = float(variance)
    print(f"\nThe numerical value of the variance is: {final_answer}")

solve_variance()