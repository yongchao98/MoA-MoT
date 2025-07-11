from fractions import Fraction

def calculate_a2_coefficients():
    """
    Calculates the coefficients of the Seeley-DeWitt a_2 term
    for a massless gauged Dirac spinor field in 4 dimensions.
    """

    # The general formula for the trace of the local a_2 coefficient for an
    # operator P = -Delta - E_P on a 4D manifold is, ignoring total derivatives:
    # tr(a_2) = tr( C1*(R_....^2 - R_..^2)*I + C2*R^2*I + C3*R*E_P + C4*E_P^2 - C5*Omega^2 )
    #
    # The standard coefficients are:
    # C1 = 1/180, C2 = 1/72, C3 = 1/6, C4 = 1/2, C5 = 1/12
    #
    # For the squared Dirac operator D^2, we have:
    # E_P = -R/4 + (1/2)sigma^munu * F_munu
    # Omega_munu = R_munu(spin) + F_munu
    #
    # We need to compute the traces of various terms over the fiber (spinor x gauge).
    # N is the dimension of the gauge representation (e.g., N for U(N)).
    # The dimension of the Dirac spinor space is 4.
    #
    # Required trace relations:
    # tr(I) = 4*N
    # tr(E_P) = -N*R
    # tr(E_P^2) = (N/4)*R^2 + 2*tr_g(F^2)
    # tr(Omega^2) = 2*N*R_....^2 + 4*tr_g(F^2)
    # where F^2 is F_munu * F^munu and tr_g is the trace over gauge indices.

    # General formula coefficients
    c_R4_term = Fraction(1, 180)
    c_R2_term = Fraction(-1, 180)
    c_Rsq_term_I = Fraction(1, 72)
    c_RE_term = Fraction(1, 6)
    c_Esq_term = Fraction(1, 2)
    c_Omega_sq_term = Fraction(1, 12)

    # Trace factors for each component
    tr_I_factor = 4
    tr_E_R_factor = -1
    tr_Esq_Rsq_factor = Fraction(1, 4)
    tr_Esq_Fsq_factor = 2
    tr_Omega_sq_R4_factor = 2
    tr_Omega_sq_Fsq_factor = 4

    # Calculate the final coefficient for each term in the a_2 expansion.
    # The result is of the form:
    # tr(a_2) = N * (C_R4 * R_....^2 + C_R2 * R_..^2 + C_Rsq * R^2) + C_F2 * tr_g(F^2)

    # Coefficient for tr_g(F^2)
    C_F2 = c_Esq_term * tr_Esq_Fsq_factor - c_Omega_sq_term * tr_Omega_sq_Fsq_factor

    # Coefficient for N * R_....^2
    C_R4 = c_R4_term * tr_I_factor - c_Omega_sq_term * tr_Omega_sq_R4_factor

    # Coefficient for N * R_..^2
    C_R2 = c_R2_term * tr_I_factor

    # Coefficient for N * R^2
    C_Rsq = (c_Rsq_term_I * tr_I_factor) + \
            (c_RE_term * tr_E_R_factor) + \
            (c_Esq_term * tr_Esq_Rsq_factor)

    print("The second non-trivial Seeley-DeWitt coefficient, a_2, determines the Yang-Mills and quadratic gravity terms.")
    print("Its fiber trace, tr(a_2), for a massless gauged Dirac spinor is given by the expression:")
    print("tr(a_2) = N * (C_R4 * R_munurhosigma^2 + C_R2 * R_munu^2 + C_Rsq * R^2) + C_F2 * tr_g(F_munu^2)")
    print("\nWhere N is the dimension of the gauge representation, and the curvatures are:")
    print("  R_munurhosigma: Riemann curvature tensor")
    print("  R_munu: Ricci curvature tensor")
    print("  R: Ricci scalar")
    print("  F_munu: Gauge field strength")
    print("\nThe calculated rational coefficients are:")
    print(f"C_R4  = {C_R4}")
    print(f"C_R2  = {C_R2}")
    print(f"C_Rsq = {C_Rsq}")
    print(f"C_F2  = {C_F2}")

    print("\nThus, the final expression for the coefficient is:")
    # Using a_2 as a shorthand for tr(a_2) as common in physics literature.
    print(f"a_2 = N * (({C_R4}) R_munurhosigma^2 + ({C_R2}) R_munu^2 + ({C_Rsq}) R^2) + ({C_F2}) tr_g(F_munu^2)")

calculate_a2_coefficients()