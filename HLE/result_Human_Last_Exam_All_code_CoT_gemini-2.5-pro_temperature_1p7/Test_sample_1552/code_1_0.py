from fractions import Fraction

def find_heat_kernel_coefficient_a2():
    """
    Calculates and prints the coefficients of the a_2 Seeley-DeWitt term
    in the heat kernel expansion for a massless gauged Dirac spinor field in 4D.

    The local density a_2(x) has the general form:
    a_2(x) = N_f * [ C_Riem^2 * (R_μνab)^2 + C_Ricci^2 * (R_μν)^2 + C_R^2 * R^2 + C_F^2 * tr(F_μν^2) ]
    where N_f is the dimension of the gauge group representation.
    """

    # The coefficients are derived from established results in quantum field theory,
    # specifically for a 4-component Dirac fermion.
    # The gravitational part is often expressed in terms of the squared Weyl tensor (C^2),
    # the Gauss-Bonnet term (E_4), and the squared scalar curvature (R^2).
    # a_2(x) / N_f = (1/120)C^2 - (1/360)E_4 + (1/72)R^2 - (1/6)tr(F^2)
    # where C^2 = R_μνab^2 - 2R_μν^2 + (1/3)R^2 and E_4 = R_μνab^2 - 4R_μν^2 + R^2.
    # We expand these terms to find the coefficients of the individual curvature invariants.

    # Coefficient for R_μνab * R^μνab
    c_riemann_sq = Fraction(1, 120) - Fraction(1, 360)

    # Coefficient for R_μν * R^μν
    c_ricci_sq = Fraction(-2, 120) + Fraction(4, 360)

    # Coefficient for R^2
    c_r_sq = Fraction(1, 360) - Fraction(1, 360) + Fraction(1, 72)
    
    # Coefficient for tr_V(F_μν * F^μν)
    # This comes from adding the contributions from the field curvature Ω and the endomorphism term X.
    # The consensus value for a Dirac fermion, after accounting for normalization conventions, is -1/6.
    c_f_sq = Fraction(-1, 6)
    
    print("The second Seeley-DeWitt coefficient density, a_2(x), for a massless gauged")
    print("Dirac spinor field in a 4-dimensional Riemannian manifold is given by an expression")
    print("proportional to the dimension of the gauge representation, N_f:")
    print("\na_2(x) = N_f * [ C_Riem^2 * (R_μνab R^μνab) + C_Ricci^2 * (R_μν R^μν) + C_R^2 * R^2 + C_F^2 * tr_V(F_μν F^μν) ]\n")
    print("The coefficients are:")
    print(f"C_Riem^2 (coefficient of R_μνab R^μνab): {c_riemann_sq}")
    print(f"C_Ricci^2 (coefficient of R_μν R^μν):    {c_ricci_sq}")
    print(f"C_R^2 (coefficient of R^2):               {c_r_sq}")
    print(f"C_F^2 (coefficient of tr_V(F_μν F^μν)):  {c_f_sq}\n")
    
    print("So, the final equation is:")
    print(f"a_2(x) = N_f * [ ({c_riemann_sq}) * R_μνab R^μνab + ({c_ricci_sq}) * R_μν R^μν + ({c_r_sq}) * R^2 + ({c_f_sq}) * tr_V(F_μν F^μν) ]")

if __name__ == '__main__':
    find_heat_kernel_coefficient_a2()