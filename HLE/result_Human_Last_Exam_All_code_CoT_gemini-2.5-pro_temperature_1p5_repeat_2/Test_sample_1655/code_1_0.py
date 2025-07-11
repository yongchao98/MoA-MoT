import sympy as sp
from sympy import Function, Matrix, Symbol

def main():
    """
    Analyzes an optical system symbolically to determine if an arbitrary vector
    beam can be generated from a linearly polarized input.
    """
    # --- Symbolic Setup ---
    # Define symbolic variables and functions for a 1D system for simplicity.
    # The same logic applies to a 2D system (x, y).
    x = Symbol('x', real=True)
    E_in = Function('E_in')(x)
    t11 = Function('t11')(x)
    t12 = Function('t12')(x)
    t21 = Function('t21')(x)
    t22 = Function('t22')(x)

    # Define P as a generic, unevaluated linear operator representing free-space propagation.
    # P(f(x)) means "propagation of the field f(x)".
    class P(Function):
        # We enforce linearity: P(a*f + b*g) = a*P(f) + b*P(g)
        @classmethod
        def eval(cls, arg):
            if arg.is_Add:
                return sp.Add(*[cls(a) for a in arg.args])
            if arg.is_Mul:
                coeffs, funcs = arg.as_coeff_mul()
                if funcs: # Ensure there is a function to operate on
                    return coeffs * cls(sp.Mul(*funcs))
            return None # Let sympy handle other cases or leave unevaluated

    # --- System Simulation ---
    print("This script symbolically analyzes an optical system.")
    print("The goal is to determine if an arbitrary vector beam can be generated from an input beam")
    print("with a fully controllable scalar field (phase and amplitude) but fixed linear polarization.\n")
    
    # 1. Input Beam
    # The input has a tailored scalar field E_in(x) and is linearly polarized along the x-axis.
    # The first free-space propagation is absorbed into the definition of E_in(x), as we have
    # full control over the field at the input of the first random medium.
    E_input_vector = Matrix([E_in, 0])
    print("Step 1: Define the input vector field (linearly polarized).")
    print(f"E_input = {E_input_vector}\n")

    # 2. First Random Medium (Matrix T)
    T = Matrix([[t11, t12], [t21, t22]])
    E_after_T = T * E_input_vector
    print("Step 2: The beam passes through the first random medium T(x).")
    print(f"Field after T = T * E_input = {E_after_T}\n")
    print("Note: Both polarization components are now determined by the *single* input function E_in(x).\n")

    # 3. Second Free-Space Propagation (Operator P)
    # The propagation operator P acts on each polarization component independently.
    E_after_P = Matrix([
        P(E_after_T[0]),
        P(E_after_T[1])
    ])
    print("Step 3: The resulting vector beam propagates through free space.")
    print("Let P[f(x)] represent the propagation of a field f(x). P is a linear operator.")
    print(f"Field after propagation P =\n{E_after_P}\n")

    # 4. Second Random Medium (Inverse Matrix T_inv)
    # The beam passes through a medium with the inverse transmission matrix.
    T_inv = T.inv()
    E_output_vector = T_inv * E_after_P

    # For clarity, let's expand the final expression manually using P's linearity
    det_T = t11*t22 - t12*t21
    E_out_x = sp.expand((1/det_T) * (t22 * E_after_P[0] - t12 * E_after_P[1]))
    E_out_y = sp.expand((1/det_T) * (-t21 * E_after_P[0] + t11 * E_after_P[1]))
    
    print("Step 4: The beam passes through the second medium, T_inv(x).")
    print("The final output field E_out = T_inv * (Field after propagation P).\n")
    
    print("--- Final Output Field Components ---")
    print("The final x-component of the electric field is:")
    sp.pprint(E_out_x, use_unicode=False)
    print("\nThe final y-component of the electric field is:")
    sp.pprint(E_out_y, use_unicode=False)
    print("\n")

    # --- Analysis and Conclusion ---
    print("--- Analysis ---")
    print("To generate an 'arbitrary vector beam', we must be able to specify two independent complex functions,")
    print("E_out_x and E_out_y, by choosing the right input.")
    print("\nHowever, the symbolic derivation above shows that both E_out_x and E_out_y are constructed")
    print("from the same two fundamental terms: P[t11(x)*E_in(x)] and P[t21(x)*E_in(x)].")
    print("These terms, in turn, both depend on the *single* controllable input function, E_in(x).")
    print("\nBecause we only have one degree of freedom (the scalar field E_in(x)) to control the input,")
    print("we cannot independently control the two output functions. If we choose E_in(x) to produce")
    print("a desired E_out_x, the resulting E_out_y is automatically determined by the system's physics")
    print("and cannot be chosen freely. Therefore, the set of achievable output beams is a restricted")
    print("subset, not the set of all arbitrary vector beams.")

if __name__ == "__main__":
    main()
<<<No>>>