import math

def solve_for_t0():
    """
    Calculates the positive value of t0 based on the solvability condition
    of the given boundary-value problem.
    """
    # Given parameters from the problem description
    alpha = 10**16
    R = math.log(100.0 / 99.0)

    # The solvability condition leads to the algebraic equation:
    # (t0**2 / 3) * (e**R - 1) = alpha
    # We solve this equation for t0.

    # Calculate the term (e^R - 1)
    eR_minus_1 = math.exp(R) - 1

    # Calculate t0 squared
    # t0^2 = 3 * alpha / (e^R - 1)
    t0_squared = 3 * alpha / eR_minus_1

    # Calculate the positive value of t0
    t0 = math.sqrt(t0_squared)

    # Output the components of the final equation and the result
    three = 3
    print(f"The final equation for t0^2 is: ({three} * alpha) / (e^R - 1)")
    print(f"where alpha = {alpha}")
    print(f"and R = {R}")
    print(f"The term (e^R - 1) evaluates to: {eR_minus_1}")
    print(f"Plugging in the numbers: t0^2 = ({three} * {alpha}) / {eR_minus_1} = {t0_squared}")
    print(f"The positive value of t0 is: {t0}")

if __name__ == "__main__":
    solve_for_t0()
