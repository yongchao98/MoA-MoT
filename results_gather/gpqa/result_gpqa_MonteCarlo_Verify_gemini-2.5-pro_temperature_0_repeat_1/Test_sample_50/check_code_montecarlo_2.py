import sympy

def check_potential_energy_formula():
    """
    Checks the correctness of the selected answer for the potential energy
    of a charge near a grounded conducting sphere.

    The function uses symbolic mathematics to compare the formula from the
    selected option with the theoretically derived correct formula.
    """
    # Define the physical quantities as symbolic variables
    # k: Coulomb's constant
    # q: charge
    # R: radius of the sphere
    # d: distance of the charge from the center of the sphere
    k, q, R, d = sympy.symbols('k q R d', real=True, positive=True)

    # The problem implies d > R (charge is outside the sphere)

    # The correct formula for the potential energy of the system, derived
    # by integrating the force from the image charge as the real charge
    # is brought from infinity.
    # U = -(1/2) * k * q^2 * R / (d^2 - R^2)
    correct_formula = -(1/2) * k * q**2 * R / (d**2 - R**2)

    # The formula from the selected answer (Option D)
    # D) U = -(1/2) * kq^2 R/(d^2 -R^2)
    answer_formula_D = -(1/2) * k * q**2 * R / (d**2 - R**2)

    # The formulas from the other options for completeness
    answer_formula_A = -k * q**2 * d / (d**2 - R**2)
    answer_formula_B = -(1/2) * k * q**2 * R**2 / (d**2 - R**2)
    answer_formula_C = -(1/2) * k * q**2 * d / (d**2 + R**2)

    # Check if the selected answer's formula is mathematically identical
    # to the theoretically correct formula.
    # We use sympy.simplify to ensure the comparison is robust.
    if sympy.simplify(correct_formula - answer_formula_D) == 0:
        # The formula is correct. Now check the derivation provided in the answer.
        # Derivation steps:
        # 1. Image charge q' = -q * (R/d)
        # 2. Image position b = R^2 / d
        # 3. Potential at q due to image: V_induced = k * q' / (d - b)
        # 4. Energy U = (1/2) * q * V_induced
        q_prime = -q * R / d
        b = R**2 / d
        V_induced = k * q_prime / (d - b)
        U_derived_in_answer = sympy.simplify((1/2) * q * V_induced)

        if sympy.simplify(U_derived_in_answer - correct_formula) == 0:
            return "Correct"
        else:
            return (f"Incorrect. The final formula in option D is correct, but the "
                    f"derivation provided in the explanation is flawed. The derived "
                    f"expression was {U_derived_in_answer}, which does not simplify "
                    f"to the correct answer.")
    else:
        # If the formula for D was wrong, identify the correct option.
        if sympy.simplify(correct_formula - answer_formula_A) == 0:
            correct_option = "A"
        elif sympy.simplify(correct_formula - answer_formula_B) == 0:
            correct_option = "B"
        elif sympy.simplify(correct_formula - answer_formula_C) == 0:
            correct_option = "C"
        else:
            correct_option = "not among the options"

        return (f"Incorrect. The selected answer is D, but the correct formula for the "
                f"potential energy is U = -(1/2) * k * q^2 * R / (d^2 - R^2). "
                f"This corresponds to option {correct_option}.")

# Execute the check and print the result
result = check_potential_energy_formula()
print(result)