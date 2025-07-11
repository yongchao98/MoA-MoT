from fractions import Fraction

def compute_scl():
    """
    This function computes the stable commutator length (scl) for the given problem.
    """
    # Step 1: Define the rotation numbers from the problem statement.
    # g is translation by 2/27, so its rotation number is 2/27.
    # h is translation by 16/27, so its rotation number is 16/27.
    rot_g = Fraction(2, 27)
    rot_h = Fraction(16, 27)

    # Step 2: Apply the formula for scl of a product in a free product group.
    # The formula is scl(g_1 * h_2) = |rot(g_1) - rot(h_2)| / 2.
    diff = rot_g - rot_h
    abs_diff = abs(diff)
    scl_value = abs_diff / 2

    # Step 3: Print the calculation step-by-step, showing each number in the equation.
    print(f"The rotation number for g is {rot_g.numerator}/{rot_g.denominator}")
    print(f"The rotation number for h is {rot_h.numerator}/{rot_h.denominator}")
    print("\nThe stable commutator length (scl) is calculated using the formula:")
    print("scl = |rotation_number(g) - rotation_number(h)| / 2\n")

    print(f"scl = |{rot_g.numerator}/{rot_g.denominator} - {rot_h.numerator}/{rot_h.denominator}| / 2")
    print(f"scl = |{diff.numerator}/{diff.denominator}| / 2")
    print(f"scl = ({abs_diff.numerator}/{abs_diff.denominator}) / 2")
    print(f"scl = {scl_value.numerator}/{scl_value.denominator}")

if __name__ == "__main__":
    compute_scl()
