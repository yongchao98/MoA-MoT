def get_critical_exponent_nu():
    """
    Provides the precise value for the critical exponent nu (ν) in a
    3-dimensional G₄ (interpreted as φ⁴) theoretical framework.
    """
    # The spatial dimension under consideration. The d=3 case is the most
    # common non-trivial example in statistical mechanics.
    dimension = 3

    # For d=3, the Ising universality class (described by φ⁴-theory) has a
    # critical exponent ν known to very high precision from field theory,
    # Monte Carlo simulations, and the conformal bootstrap.
    # We use a state-of-the-art value here.
    nu = 0.629971

    # The correlation length ξ is a measure of the typical size of correlated
    # fluctuations. Near the critical temperature T_c, it diverges according
    # to a power law.
    # The request to "output each number in the final equation" is interpreted
    # as showing the scaling relation with its components.
    
    print(f"Assuming the 'G₄-theoretical framework' refers to the φ⁴-theory (Ising universality class) in d = {dimension} spatial dimensions:")
    print(f"The scaling of the correlation length ξ near the critical temperature T_c is given by the relation:")
    print(f"ξ ∝ |T - T_c|^(-ν)")
    print("\nIn this equation:")
    print(f"- The exponent ν is the critical exponent for the correlation length.")
    print(f"- Its precise value for d = {dimension} is ν = {nu}")

if __name__ == "__main__":
    get_critical_exponent_nu()
