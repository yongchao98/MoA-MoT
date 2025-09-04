import sympy

def check_potential_energy_formula():
    """
    Checks the correctness of the formula for the potential energy of a charge
    near a grounded conducting sphere.

    The function performs two main checks:
    1. Derives the correct formula from the method of images using symbolic math.
    2. Verifies that the chosen answer satisfies critical physical limits.
    """
    try:
        # 1. Define symbolic variables for the physical quantities.
        # We assume all are positive real numbers.
        k, q, R, d = sympy.symbols('k q R d', positive=True, real=True)

        # 2. Define the candidate answer to be checked.
        # The provided answer is A: U = - (1/2) * k * q^2 * R / (d^2 - R^2)
        llm_answer_formula = -(1/2) * k * q**2 * R / (d**2 - R**2)

        # 3. Derive the correct formula from the method of images.
        # Image charge magnitude: q' = -q * (R/d)
        q_image = -q * R / d
        # Image charge position from center: b = R^2 / d
        b_image = R**2 / d
        # Distance between the real charge q and the image charge q'
        distance_between_charges = d - b_image
        # The potential at the location of q due to the image charge q'
        V_image = k * q_image / distance_between_charges
        # The net potential energy of the system is (1/2) * q * V_image
        derived_formula = sympy.simplify((1/2) * q * V_image)

        # 4. Compare the LLM's answer with the derived formula.
        # The difference should be zero if they are symbolically identical.
        if sympy.simplify(llm_answer_formula - derived_formula) != 0:
            return (f"Incorrect. The provided answer's formula is {llm_answer_formula}, "
                    f"but the correct formula derived from the method of images is {derived_formula}.")

        # 5. Check if the formula satisfies physical constraints (as a sanity check).
        # Constraint 1: As the charge approaches the sphere (d -> R+), the potential
        # energy should go to negative infinity due to the infinite attractive force.
        limit_d_to_R = sympy.limit(llm_answer_formula, d, R, dir='+')
        if limit_d_to_R != -sympy.oo:
            return (f"Incorrect. The formula does not satisfy the physical constraint "
                    f"that U -> -infinity as d -> R. The calculated limit is {limit_d_to_R}.")

        # Constraint 2: As the charge moves very far away (d -> infinity), the potential
        # energy should approach zero.
        limit_d_to_inf = sympy.limit(llm_answer_formula, d, sympy.oo)
        if limit_d_to_inf != 0:
            return (f"Incorrect. The formula does not satisfy the physical constraint "
                    f"that U -> 0 as d -> infinity. The calculated limit is {limit_d_to_inf}.")
        
        # Constraint 3: For large d, the sphere acts like an induced dipole. The potential
        # energy of a charge-dipole interaction falls off as 1/d^2.
        series_expansion = sympy.series(llm_answer_formula, d, sympy.oo, n=3)
        # The leading term in the series for large d should be proportional to 1/d^2.
        # The `getn(2)` method gets the coefficient of the d**(-2) term.
        # If it's zero, the fall-off is faster. If it's non-zero, the fall-off is 1/d^2.
        # We check that the order of the leading term is d**-2.
        if sympy.O(1/d**2, (d, sympy.oo)) not in series_expansion.args:
             return (f"Incorrect. The potential energy should fall off as 1/d^2 for large d. "
                     f"The series expansion is {series_expansion}, which indicates a different behavior.")


        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_potential_energy_formula()
print(result)