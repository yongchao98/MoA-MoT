import sympy

def calculate_yukawa_ratio():
    """
    Calculates the ratio R = delta_Zx / (delta_Zg + delta_Zm_x) for a Yukawa theory.

    The calculation is based on the one-loop counter-terms in the MS-bar scheme.
    We represent the counter-terms as numerical coefficients times a common factor
    C = g^2 / (32 * pi^2 * epsilon), which will cancel in the final ratio.
    """

    # 1. Fermion wave-function counter-term coefficient (delta_Zx)
    # From the fermion self-energy diagram, the coefficient of the divergent part
    # corresponding to the wave-function renormalization is -1.
    coeff_zx = -1
    print(f"The coefficient for the fermion field counter-term, delta_Zx, is proportional to: {coeff_zx}")


    # 2. Fermion mass counter-term coefficient (delta_Zm_x)
    # From the fermion self-energy, the coefficient for the mass counter-term is 3.
    coeff_zmx = 3
    print(f"The coefficient for the fermion mass counter-term, delta_Zm_x, is proportional to: {coeff_zmx}")

    # 3. Yukawa coupling counter-term coefficient (delta_Zg)
    # This is derived from the relation: delta_Zg = delta_1 - delta_Zx - (1/2)*delta_Zphi
    # where delta_1 is the vertex counter-term and delta_Zphi is the scalar field counter-term.

    # For a scalar Yukawa theory, the vertex counter-term is equal to the fermion wave-function counter-term.
    coeff_d1 = coeff_zx
    print(f"The coefficient for the vertex part, delta_1, is equal to that of delta_Zx: {coeff_d1}")

    # The problem states that the scalar field renormalization constant is zero.
    coeff_zphi = 0
    print(f"The problem states delta_Zphi = 0, so its coefficient is: {coeff_zphi}")

    # Now, calculate the coefficient for delta_Zg
    coeff_zg = coeff_d1 - coeff_zx - 0.5 * coeff_zphi
    print(f"The coefficient for the Yukawa coupling counter-term, delta_Zg, is calculated as delta_1 - delta_Zx - 0.5*delta_Zphi = {coeff_d1} - ({coeff_zx}) - 0.5*({coeff_zphi}) = {coeff_zg}")

    # 4. Calculate the final ratio R
    # The common factor C cancels out, so we only need the coefficients.
    numerator = coeff_zx
    denominator = coeff_zg + coeff_zmx

    R = sympy.Rational(numerator, denominator)

    print("\nCalculating the ratio R = delta_Zx / (delta_Zg + delta_Zm_x)")
    print(f"R = {coeff_zx} / ({coeff_zg} + {coeff_zmx})")
    print(f"R = {numerator} / {denominator}")
    print(f"The final ratio R is: {R}")

if __name__ == "__main__":
    calculate_yukawa_ratio()
