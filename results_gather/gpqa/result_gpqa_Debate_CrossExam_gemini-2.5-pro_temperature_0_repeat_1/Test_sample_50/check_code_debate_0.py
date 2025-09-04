import sympy as sp

def check_potential_energy_formula():
    """
    Checks the correctness of the formula for the potential energy of a charge
    near a grounded conducting sphere.

    The problem:
    A charge q is placed a distance d away from the center of a grounded
    conducting sphere of radius R. Calculate the net potential energy of this system.

    The proposed answer (A):
    U = - (1/2) * k * q^2 * R / (d^2 - R^2)
    """
    # Define symbolic variables for the physical quantities.
    # k is Coulomb's constant. All variables are positive real numbers.
    k, q, R, d = sp.symbols('k q R d', real=True, positive=True)

    # The formula from the given answer A
    answer_formula_A = -(sp.S(1)/2) * k * q**2 * R / (d**2 - R**2)

    # --- Step 1: Derive the correct formula using the method of images ---
    # The potential energy U is given by U = (1/2) * q * V_induced, where
    # V_induced is the potential at the location of charge q due to the
    # induced charge on the sphere.
    # The method of images models the induced charge with an image charge q_prime
    # at a position b from the center.

    # Magnitude and position of the image charge
    q_prime = -q * R / d
    b = R**2 / d

    # The distance between the real charge q (at d) and the image charge q_prime (at b)
    distance_between_charges = d - b

    # The potential at the location of q due to the image charge (V_induced)
    V_induced = k * q_prime / distance_between_charges

    # Substitute expressions for q_prime and b to get V_induced in terms of q, R, d
    V_induced = V_induced.subs(b, R**2 / d)
    
    # Calculate the potential energy U
    derived_U = sp.S(1)/2 * q * V_induced
    
    # Simplify the derived expression
    derived_U_simplified = sp.simplify(derived_U)

    # --- Step 2: Compare the derived formula with the given answer ---
    # The difference should be zero if they are equivalent.
    if sp.simplify(derived_U_simplified - answer_formula_A) != 0:
        return (f"Incorrect. The provided formula does not match the one derived from the method of images.\n"
                f"Provided formula: {answer_formula_A}\n"
                f"Correctly derived formula: {derived_U_simplified}")

    # --- Step 3: Check physical constraints and limiting cases ---
    # Constraint: The charge must be outside the sphere, so d > R.
    # The denominator (d**2 - R**2) is positive. Since all other terms are positive,
    # the potential energy U is negative. This is physically correct, as the force
    # between the point charge and the grounded sphere is attractive.

    # Limiting case 1: The charge is infinitely far away (d -> oo).
    # The potential energy should approach 0.
    limit_d_inf = sp.limit(answer_formula_A, d, sp.oo)
    if limit_d_inf != 0:
        return (f"Incorrect. The formula fails the limiting case for d -> infinity. "
                f"Expected U -> 0, but the formula gives {limit_d_inf}.")

    # Limiting case 2: The charge approaches the surface of the sphere (d -> R+).
    # The attractive force becomes infinite, so the potential energy should go to -oo.
    limit_d_R = sp.limit(answer_formula_A, d, R, dir='+')
    if limit_d_R != -sp.oo:
        return (f"Incorrect. The formula fails the limiting case for d -> R. "
                f"Expected U -> -oo, but the formula gives {limit_d_R}.")

    # --- Step 4: Check dimensional consistency ---
    # [k*q^2] has units of [Energy * Length].
    # [R] has units of [Length].
    # [d^2 - R^2] has units of [Length^2].
    # So, the units of the formula are [Energy * Length * Length] / [Length^2] = [Energy].
    # The formula is dimensionally correct.

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_potential_energy_formula()
print(result)