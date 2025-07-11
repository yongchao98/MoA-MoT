import sympy

def find_heat_kernel_coefficient():
    """
    Calculates the second coefficient in the heat kernel expansion for a
    massless gauged Dirac spinor field.
    The calculation provides the numerical coefficient multiplying the scalar
    curvature R in the integrand of the a_2 Seeley-DeWitt coefficient.
    This coefficient is computed per dimension of the gauge group representation.
    """

    # The coefficient of R in the endomorphism E from the Lichnerowicz formula.
    coeff_E = sympy.Rational(1, 4)

    # The universal coefficient of R in the general formula for a_2.
    coeff_a2_formula = sympy.Rational(1, 6)

    # The dimension of the Dirac spinor space in 4D.
    dim_spinor = 4

    print("Step 1: Determine the coefficient of the scalar curvature R inside the trace.")
    print(f"From the Lichnerowicz formula for D^2, the endomorphism E contains a term ({coeff_E}) * R.")
    print(f"The general formula for the a_2 coefficient has a subtraction term (-{coeff_a2_formula}) * R.")
    
    coeff_inside_trace = coeff_E - coeff_a2_formula
    
    print(f"The combined coefficient for R inside the trace is: {coeff_E} - {coeff_a2_formula} = {coeff_inside_trace}")
    print("\nStep 2: Multiply by the dimension of the spinor space to perform the trace.")
    print(f"The dimension of a Dirac spinor in 4D is {dim_spinor}.")

    final_coefficient = coeff_inside_trace * dim_spinor

    print("\nFinal equation for the coefficient (per dimension of gauge representation):")
    # Print out each number and the operations in the final equation.
    print(f"({coeff_E} - {coeff_a2_formula}) * {dim_spinor} = {coeff_inside_trace} * {dim_spinor} = {final_coefficient}")

find_heat_kernel_coefficient()