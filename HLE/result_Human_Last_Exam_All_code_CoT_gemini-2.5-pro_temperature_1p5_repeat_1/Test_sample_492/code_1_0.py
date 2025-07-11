import math

def solve_for_critical_time():
    """
    Calculates the critical time 'c' for the emergence of the giant component
    in the described graph model.

    The derivation shows that at the critical time c, the following equation holds:
    c^2 / 3 = 1
    This function solves this equation for c and prints the steps.
    """

    # The final equation derived from the model is c^2 / 3 = 1.
    # We will represent the components of this equation.
    # The variable c is what we are solving for.
    numerator_term = "c^2"
    denominator = 3
    rhs = 1 # Right Hand Side

    print("The final equation relating the critical time 'c' to the system parameters simplifies to:")
    print(f"{numerator_term} / {denominator} = {rhs}")

    # To solve for c, we first solve for c^2.
    # c^2 = 3 * 1
    c_squared = float(denominator * rhs)
    print("\nSolving for c^2:")
    print(f"c^2 = {int(c_squared)}")

    # Now, we solve for c by taking the square root.
    c = math.sqrt(c_squared)
    print("\nSolving for c:")
    print(f"c = sqrt({int(c_squared)})")

    print(f"\nThe exact value of c is sqrt(3).")
    print(f"The numerical value of c is approximately: {c}")

if __name__ == "__main__":
    solve_for_critical_time()
