import sympy

def solve_quartic_euler_characteristic():
    """
    Calculates the orbifold Euler characteristic of the moduli stack of smooth
    plane quartic curves.
    """
    # The genus g of a smooth plane quartic curve (degree d=4) is (d-1)(d-2)/2.
    g = 3

    # The orbifold Euler characteristic of the moduli stack of genus g curves,
    # chi_orb(M_g), is given by the Harer-Zagier theorem as zeta(1-2g) for g>=2.
    # For g=3, this is zeta(-5).
    chi_orb_Mg = sympy.zeta(1 - 2 * g)

    # The orbifold Euler characteristic of the moduli stack of hyperelliptic
    # genus g curves, chi_orb(H_g), is given by a theorem of Bergström as
    # zeta(1-2g) for g>=2. For g=3, this is also zeta(-5).
    chi_orb_Hg = sympy.zeta(1 - 2 * g)

    # The moduli stack of interest, of smooth plane quartics, is the stack of
    # non-hyperelliptic genus 3 curves (M_3^nh). Its Euler characteristic is
    # the difference between the total space M_3 and the hyperelliptic locus H_3.
    # chi(M_3^nh) = chi(M_3) - chi(H_3)
    chi_orb_result = chi_orb_Mg - chi_orb_Hg

    # Print the step-by-step calculation
    print("Step 1: Determine the genus of a smooth plane quartic curve.")
    print(f"A plane curve of degree d=4 has genus g = (4-1)(4-2)/2 = {g}.\n")

    print("Step 2: State the formula for the Euler characteristic of the moduli space of genus 3 curves, M_3.")
    print(f"chi_orb(M_3) = zeta(1 - 2*g) = zeta({1 - 2*g})")
    print(f"zeta({1 - 2*g}) = {chi_orb_Mg}\n")

    print("Step 3: State the formula for the Euler characteristic of the hyperelliptic locus, H_3.")
    print(f"chi_orb(H_3) = zeta(1 - 2*g) = zeta({1 - 2*g})")
    print(f"zeta({1 - 2*g}) = {chi_orb_Hg}\n")

    print("Step 4: The final result is the difference between these two values.")
    print("The final equation is: chi([U/G]) = chi_orb(M_3) - chi_orb(H_3)")
    print(f"Result = {chi_orb_Mg} - ({chi_orb_Hg}) = {chi_orb_result}")


solve_quartic_euler_characteristic()