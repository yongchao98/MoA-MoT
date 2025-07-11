import math

def solve_exponent_order():
    """
    This function explains and determines the order in the coupling constant 'u'
    at which the critical exponent nu first receives a non-vanishing contribution
    in phi^4 theory.
    """

    # In the renormalization group analysis of φ⁴ theory, the critical exponent ν
    # is related to the mass anomalous dimension, γ_m(u), at the Wilson-Fisher fixed point.
    # The relationship is ν = 1 / (2 - γ_m(u)).

    # The mean-field value of ν is obtained at the Gaussian fixed point, where the coupling u = 0.
    # At this point, γ_m(0) = 0.
    u_gaussian = 0
    gamma_m_at_gaussian = 0
    nu_mf_numerator = 1
    nu_mf_denominator = 2
    nu_mean_field = nu_mf_numerator / nu_mf_denominator

    print("The critical exponent ν is calculated using the formula:")
    print("ν = 1 / (2 - γ_m(u))")
    print("where u is the coupling constant and γ_m(u) is the mass anomalous dimension.")
    print("\nIn mean-field theory (at the Gaussian fixed point where u = 0), γ_m(0) = 0. This gives:")
    print(f"ν_0 = {nu_mf_numerator} / ({nu_mf_denominator} - {gamma_m_at_gaussian}) = {nu_mean_field}")

    # To find contributions beyond mean-field theory, γ_m(u) is calculated perturbatively
    # as a power series in the coupling constant u.
    # The expansion takes the form: γ_m(u) = c₁·u¹ + c₂·u² + ...

    # The first non-vanishing contribution to ν (the correction to the mean-field value)
    # comes from the lowest-order term in this series with a non-zero coefficient.
    # Standard one-loop calculations in the ε-expansion show that the coefficient c₁ is non-zero.
    order_of_first_contribution = 1

    print("\nCorrections to the mean-field value are found by calculating γ_m(u) at the non-trivial fixed point.")
    print("Perturbative calculations expand γ_m(u) in powers of u:")
    print("γ_m(u) = c₁·u¹ + c₂·u² + ...")
    print("One-loop calculations show that the first coefficient, c₁, is non-zero.")
    print("\nTherefore, the initial non-vanishing contribution to ν appears at the first order in the coupling constant u.")
    print("\nFinal equation for the order (n):")
    print(f"Contribution ∝ uⁿ, where n = {order_of_first_contribution}")


solve_exponent_order()
<<<1>>>