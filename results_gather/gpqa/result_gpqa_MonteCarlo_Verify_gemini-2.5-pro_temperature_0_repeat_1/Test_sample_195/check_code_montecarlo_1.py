import sympy

def check_answer():
    """
    Verifies the correct formula for the maximum speed of a 1D relativistic
    harmonic oscillator using symbolic mathematics.
    """
    try:
        # Define symbolic variables. Assume all are positive real numbers.
        m, k, A, c = sympy.symbols('m k A c', positive=True, real=True)

        # --- Step 1: Derive the correct expression from first principles ---

        # Energy at maximum amplitude (x=A, v=0) is purely potential
        E_potential_max = sympy.Rational(1, 2) * k * A**2

        # By conservation of energy, this equals the max kinetic energy:
        # E_potential_max = (gamma_max - 1) * m * c**2

        # Solve for gamma_max
        gamma_max_expr = 1 + E_potential_max / (m * c**2)
        # This simplifies to: 1 + (k * A**2) / (2 * m * c**2)

        # The relationship between v_max and gamma_max is v_max = c * sqrt(1 - 1/gamma_max^2)
        v_max_derived = c * sympy.sqrt(1 - 1 / gamma_max_expr**2)

        # --- Step 2: Define the expression for the given answer (Option C) ---
        v_max_C = c * sympy.sqrt(1 - 1 / (1 + (k * A**2) / (2 * m * c**2))**2)

        # --- Step 3: Verify if the derived expression is equivalent to Option C ---
        # We simplify their difference. If it's zero, they are mathematically identical.
        if sympy.simplify(v_max_derived - v_max_C) != 0:
            return "Incorrect. The formula for Option C does not match the one derived from first principles."

        # --- Step 4: Perform sanity checks on other options ---

        # Option A: v_max = c*sqrt(1 + 1/(1 - kA^2/(2mc^2)))
        # For the expression to be real, 1 - kA^2/(2mc^2) > 0.
        # In this case, the term under the square root is > 1, which implies v_max > c.
        # This is physically impossible for a massive particle.
        
        # Option B: v_max = c*sqrt(1 + 1/(1 - kA^2/(2m))^2)
        # The term k*A^2/m has units of velocity squared, while 1 is dimensionless.
        # Subtracting a value with units from a dimensionless number is a dimensional inconsistency.

        # Option D: v_max = sqrt(kA^2/m)
        # This is the classical (non-relativistic) result. We check if Option C reduces to this.
        # We find the Taylor series expansion of Option C for c -> infinity.
        classical_limit_of_C = sympy.series(v_max_C, c, sympy.oo, 3).lseries(c).next()[0]
        classical_formula_D = sympy.sqrt(k * A**2 / m)
        
        if sympy.simplify(classical_limit_of_C - classical_formula_D) != 0:
            return f"Incorrect. The non-relativistic limit of Option C is wrong. Expected {classical_formula_D}, but got {classical_limit_of_C}."

        # --- Conclusion ---
        # The formula for C is identical to the one derived from first principles,
        # and it has the correct non-relativistic limit. The other options are flawed.
        return "Correct"

    except ImportError:
        return "Verification failed: The 'sympy' library is required. Please install it using 'pip install sympy'."
    except Exception as e:
        return f"An unexpected error occurred during verification: {e}"

# Run the check
result = check_answer()
if result == "Correct":
    print("Correct")
else:
    print(f"Incorrect. Reason: {result}")
