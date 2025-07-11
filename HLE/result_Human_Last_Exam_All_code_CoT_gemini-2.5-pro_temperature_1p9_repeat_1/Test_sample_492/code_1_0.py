import math

def solve_for_c():
    """
    This function outlines the derivation for the critical time 'c' and prints the solution.

    The derivation finds that the average degree 'd' of the graph at time 'c' is d = c^2 / 3.
    The giant component emerges when this average degree d equals 1.
    This function solves the resulting equation: c^2 / 3 = 1.
    """

    # The critical equation is d = 1, where d = c^2 / 3.
    # We want to print each number in the final equation: c^2 / 3 = 1
    # The numbers are the exponent '2', the denominator '3', and the result '1'.
    base = 'c'
    power = 2
    denominator = 3
    result = 1

    # Print the explanation and the setup of the equation
    print("Let c be the time of the emergence of the giant component.")
    print("The analysis shows that the average vertex degree 'd' in the graph at time c is given by the function: d(c) = c^2 / 3.")
    print("\nA giant component emerges when the average degree is equal to 1.")
    print("Therefore, we need to solve the following equation for c:")

    # Print the equation with its numerical components
    print(f"\n({base}^{power}) / {denominator} = {result}")

    # Solve the equation step-by-step
    multiplied_result = result * denominator
    print(f"\nFirst, multiply both sides by {denominator}:")
    print(f"{base}^{power} = {multiplied_result}")

    print(f"\nNext, take the square root of both sides.")
    print("Since time c must be positive, we take the positive root:")
    print(f"{base} = sqrt({multiplied_result})")

    # Calculate and print the final exact value
    final_value_str = f"sqrt({multiplied_result})"
    final_value_numeric = math.sqrt(multiplied_result)

    print(f"\nThe exact value for c is {final_value_str}.")
    print(f"Numerically, this is approximately {final_value_numeric:.4f}.")

if __name__ == "__main__":
    solve_for_c()
