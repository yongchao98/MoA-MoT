import sympy

def mori_tanaka_model_expression():
    """
    This function generates and prints the symbolic expression for the
    effective elastic moduli of a composite material using the Mori-Tanaka model.
    """
    # Set up pretty printing for a more readable mathematical output in capable terminals.
    # In environments where this is not supported, standard printing will be used.
    try:
        sympy.init_printing(use_unicode=True)
    except ImportError:
        pass # Continue without pretty printing if not available

    # Step 1: Define the symbolic variables for the model.

    # Tensors are defined as non-commutative symbols, as the order of multiplication
    # (tensor dot product) is important.
    C, Cm, Cf, I, A = sympy.symbols('C Cm Cf I A', commutative=False)

    # Volume fractions (Vf, Vm) are scalars, so they are commutative.
    Vf, Vm = sympy.symbols('Vf Vm')

    # Step 2: Construct the Mori-Tanaka equation for the effective stiffness tensor C.
    # The model defines C as a correction to the matrix stiffness Cm.
    # The term (Vm*I + Vf*A) is inverted, which is represented by raising it to the power of -1.
    # Sympy handles the non-commutative multiplication order correctly.
    
    # We explicitly define the term to be inverted for clarity.
    term_to_invert = (Vm * I + Vf * A)
    
    # We construct the final expression for C.
    expression_for_C = Cm + Vf * (Cf - Cm) * A * term_to_invert**(-1)

    # Create a symbolic equation object for a clean representation.
    final_equation = sympy.Eq(C, expression_for_C)

    # Step 3: Print the final results to the user.
    print("The expression for the effective average elastic moduli, C, from the Mori-Tanaka model is:")
    
    # Use pretty print if available, otherwise use standard print.
    try:
        print(sympy.pretty(final_equation))
    except NameError:
        print(final_equation)


    print("\nIn this equation, each symbol represents the following quantity:")
    print(f"  C: The effective average elastic moduli tensor of the composite.")
    print(f" Cm: The fourth-order elasticity tensor of the polymer matrix.")
    print(f" Cf: The fourth-order elasticity tensor of the fiber.")
    print(f" Vf: The volume fraction of the fibers.")
    print(f" Vm: The volume fraction of the matrix (where Vm = 1 - Vf).")
    print(f"  I: The fourth-order identity tensor.")
    print(f"  A: The Eshelby strain-concentration tensor.")

if __name__ == '__main__':
    mori_tanaka_model_expression()