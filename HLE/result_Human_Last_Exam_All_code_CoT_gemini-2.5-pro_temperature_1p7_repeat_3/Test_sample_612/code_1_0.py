import sympy

def generate_mori_tanaka_expression():
    """
    This script uses the SymPy library to symbolically represent and construct
    the expression for the effective average elastic moduli (C) of a fiber-reinforced
    composite based on the Mori–Tanaka model.
    """
    # --- Step 1: Define the symbolic quantities ---
    # We treat the fourth-order tensors as non-commutative symbolic objects,
    # as tensor multiplication (contraction) is not generally commutative.
    I = sympy.Symbol('I', commutative=False)
    Cf = sympy.Symbol('C_f', commutative=False)
    Cm = sympy.Symbol('C_m', commutative=False)
    A = sympy.Symbol('A', commutative=False)

    # Vf is the scalar volume fraction of the fibers.
    Vf = sympy.Symbol('V_f')
    
    # The target variable is C, the effective elastic moduli.
    C_symbol = sympy.Symbol('C', commutative=False)

    # The volume fraction of the matrix, Vm, is derived from Vf.
    Vm = 1 - Vf

    # --- Step 2: Construct the Mori-Tanaka expression ---
    # The model provides an expression for C as the matrix stiffness plus a
    # correction term that accounts for the reinforcing fibers.
    # The formula is: C = Cm + Vf * (Cf - Cm) * A * [Vm*I + Vf*A]^(-1)
    # In this symbolic representation:
    #   - '*' represents tensor contraction (e.g., a double dot product).
    #   - 'X**-1' represents the inverse of tensor X.

    # We build the expression step by step for clarity.
    stiffness_difference = Cf - Cm
    denominator_term = Vm * I + Vf * A
    
    # The complete expression for C.
    mori_tanaka_C = Cm + Vf * stiffness_difference * A * (denominator_term**-1)

    # --- Step 3: Display the final equation ---
    # We create an equation object for a clear and formatted output.
    final_equation = sympy.Eq(C_symbol, mori_tanaka_C)

    print("The expression for the effective average elastic moduli C, according to the Mori–Tanaka model, is:")
    
    # sympy.pretty_print provides a nicely formatted, human-readable mathematical output.
    # This prints each symbol in the final equation as requested.
    sympy.pretty_print(final_equation)

# Execute the function to display the result.
generate_mori_tanaka_expression()