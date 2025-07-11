import math

def solve_asymptotic_ode():
    """
    Calculates the parameters of the asymptotic solution for the given ODE.
    The solution is assumed to be of the form y(x) = A * (x_0 - x)^alpha.
    """
    
    # As derived from the dominant balance analysis, the exponent alpha is 1/3.
    alpha_numerator = 1
    alpha_denominator = 3
    alpha = alpha_numerator / alpha_denominator

    # As derived from the coefficient balance, A^3 = -30.
    A_cubed = -30
    
    # Calculate A
    A = math.copysign(1.0, A_cubed) * (abs(A_cubed))**alpha

    # Round A to two decimal places as requested
    A_rounded = round(A, 2)
    
    # Print the parameters
    print(f"The analysis leads to the following parameters for the asymptotic solution y(x) = A * (x_0 - x)^alpha:")
    print(f"alpha = {alpha_numerator}/{alpha_denominator}")
    print(f"A = (-30)^(1/3) which is approximately {A:.4f}")
    print(f"Rounding A to two decimal places gives {A_rounded}")
    print("\n")
    
    # Construct and print the final expression for y(x).
    # The term x_0 is an integration constant determined by the initial conditions.
    # The output format follows the requirement to show each number in the equation.
    print("The final analytical expression that approximates the solution in the large x regime (approaching the singularity) is:")
    
    # Creating the equation string with the required numbers
    final_equation = f"y(x) = {A_rounded} * (x_0 - x)^({alpha_numerator}/{alpha_denominator})"
    
    print(final_equation)

# Execute the function to find and print the solution
solve_asymptotic_ode()