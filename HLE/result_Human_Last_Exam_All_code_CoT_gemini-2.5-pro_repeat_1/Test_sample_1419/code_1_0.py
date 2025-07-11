import math

def solve_fixed_point_coupling():
    """
    This function derives and prints the leading order expression for the
    Wilson-Fisher fixed point coupling u* in phi^4 theory near four dimensions.
    """
    print("Derivation of the Wilson-Fisher Fixed Point Coupling u*")
    print("=" * 60)

    # Step 1: State the Beta Function
    print("\nStep 1: The Beta Function")
    print("In quantum field theory, the renormalization group (RG) beta function, β(u), describes the running of the coupling constant u with the energy scale.")
    print("For the ϕ^4 theory in d = 4 - ε dimensions, the one-loop beta function for the dimensionless coupling u is:")
    print("\n  β(u) = -ε * u + (3 / (16 * π²)) * u² + O(u³)\n")
    print("Here, ε = 4 - d is a small positive parameter.")
    print("-" * 60)

    # Step 2: The Fixed Point Condition
    print("\nStep 2: The Fixed Point Condition")
    print("A fixed point, u*, is a value of the coupling where the theory becomes scale-invariant.")
    print("This occurs when the beta function is zero:")
    print("\n  β(u*) = 0\n")
    print("-" * 60)

    # Step 3: Solving for the Fixed Point
    print("\nStep 3: Solving the Equation")
    print("We set the beta function to zero to find the fixed points:")
    print("\n  -ε * u* + (3 / (16 * π²)) * (u*)² = 0\n")
    print("This equation can be factored:")
    print("\n  u* * [-ε + (3 / (16 * π²)) * u*] = 0\n")
    print("This gives two solutions:")
    print("1. u* = 0. This is the trivial 'Gaussian' fixed point, corresponding to a non-interacting theory.")
    print("2. -ε + (3 / (16 * π²)) * u* = 0. This gives the non-trivial 'Wilson-Fisher' fixed point.")
    print("-" * 60)

    # Step 4: The Final Expression
    print("\nStep 4: The Leading Order Expression for u*")
    print("We solve the second case for u*:")
    print("\n  (3 / (16 * π²)) * u* = ε\n")
    print("Rearranging the terms to solve for u*, we get the final expression:")

    # Define the components of the expression
    coefficient_numerator = 16
    coefficient_denominator = 3
    pi_symbol = "π"
    epsilon_symbol = "ε"

    print("\n" + "="*25 + " FINAL RESULT " + "="*25)
    print(f"      u* = ( {coefficient_numerator} * {pi_symbol}² / {coefficient_denominator} ) * {epsilon_symbol}")
    print("=" * 64 + "\n")

if __name__ == '__main__':
    solve_fixed_point_coupling()