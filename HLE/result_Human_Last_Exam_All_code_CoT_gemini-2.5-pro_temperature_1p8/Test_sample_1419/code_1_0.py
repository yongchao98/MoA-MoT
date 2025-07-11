import math

def display_fixed_point_coupling():
    """
    This function explains and displays the leading order expression for the
    Wilson-Fisher fixed point coupling (u*) in phi^4 theory near d=4 dimensions.
    """
    print("This script derives the leading order expression for the fixed point coupling u* in terms of epsilon = 4-d.")
    print("-" * 70)

    # 1. State the Beta function
    print("\nStep 1: The one-loop beta function for the dimensionless coupling 'u' is:")
    print("β(u) = -epsilon * u + (3 / (16 * π^2)) * u^2")

    # 2. State the fixed point condition
    print("\nStep 2: A fixed point u* is found by setting the beta function to zero:")
    print("β(u*) = 0  =>  -epsilon * u* + (3 / (16 * π^2)) * (u*)^2 = 0")

    # 3. Solve for the non-trivial fixed point
    print("\nStep 3: We solve for the non-trivial (Wilson-Fisher) fixed point, where u* is not zero:")
    print("Factoring out u*: u* * [-epsilon + (3 / (16 * π^2)) * u*] = 0")
    print("The non-trivial solution comes from the second term:")
    print("(3 / (16 * π^2)) * u* = epsilon")

    # 4. Present the final expression
    print("\nStep 4: The final leading order expression for the fixed point coupling is:")

    numerator = 16
    denominator = 3
    pi_squared_symbol = "π^2"

    # Using an f-string to explicitly show the numbers in the equation
    final_equation = f"u* = ({numerator} * {pi_squared_symbol} / {denominator}) * epsilon"
    print(final_equation)

    # 5. Calculate and display the numerical value of the prefactor
    numerical_prefactor = (numerator * math.pi**2) / denominator
    print(f"\nThe numerical value of the prefactor ({numerator} * π^2 / {denominator}) is approximately {numerical_prefactor:.4f}.")
    print(f"So, numerically: u* ≈ {numerical_prefactor:.4f} * epsilon")
    print("-" * 70)

# Run the function to display the derivation
display_fixed_point_coupling()