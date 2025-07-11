import math

def solve_probability():
    """
    This function calculates the probability that the given series converges.
    It iterates through all possible values for X, Y, and Z, checks the convergence
    condition, and then computes the probability.
    """
    # Define the ranges for X, Y, and Z based on the problem statement.
    # X is a non-zero integer in [-9, 9].
    X_values = list(range(-9, 0)) + list(range(1, 10))
    # Y and Z are integers in [0, 9].
    Y_values = range(10)
    Z_values = range(10)

    # Calculate the total number of possible combinations.
    total_cases = len(X_values) * len(Y_values) * len(Z_values)

    # Initialize a counter for the number of cases where the series converges.
    convergent_cases = 0

    # Iterate through all possible combinations of X, Y, and Z.
    for X in X_values:
        for Y in Y_values:
            for Z in Z_values:
                # Calculate the variable u. The number XYZ is interpreted as 100*X + 10*Y + Z.
                # u = (100*X + 10*Y + Z)/100 + Z/10
                u = X + (10 * Y + 11 * Z) / 100.0

                # Calculate the common ratio r of the geometric series.
                r = 20 * u**2 + 24 * u

                # The series converges if the absolute value of the common ratio is less than 1.
                if abs(r) < 1:
                    convergent_cases += 1
    
    # The problem asks to output the numbers in the final equation.
    # The final equation for the probability is P = convergent_cases / total_cases.
    print("The final equation for the probability is:")
    print(f"{convergent_cases} / {total_cases}")

solve_probability()