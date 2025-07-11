import sympy

def calculate_proportionality_factor(d_val, k_val):
    """
    Calculates the proportionality factor for the given gamma matrix identity.

    The identity is:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}

    The proportionality factor C is given by:
    C(d, k) = (d - 2k)^2 - (-1)^k * (d - 2k)

    Args:
        d_val (int): The number of spacetime dimensions.
        k_val (int): The number of antisymmetrized gamma matrices in the middle term.

    Returns:
        The integer value of the proportionality factor.
    """
    if k_val < 0 or d_val < 0:
        raise ValueError("d and k must be non-negative integers.")
    if k_val > d_val:
        # In d dimensions, the antisymmetrized product of more than d gamma matrices is zero.
        # So the expression is 0 = C * 0, which means C can be anything.
        # However, typically the factor is considered in a context where gamma_k is not zero.
        # If we take the formula as primary, it still gives a value.
        pass

    d = sympy.Symbol('d')
    k = sympy.Symbol('k')

    factor_expr = (d - 2*k)**2 - (-1)**k * (d - 2*k)
    
    # Substitute the numerical values for d and k
    factor_val = factor_expr.subs({d: d_val, k: k_val})

    # Print the equation with the calculated factor
    print(f"For d = {d_val} and k = {k_val}, the equation is:")
    print(f"γ_μν γ_{'μ₁...μₖ' if k_val > 1 else 'μ₁'} γ^μν = {factor_val} γ_{'μ₁...μₖ' if k_val > 1 else 'μ₁'}")
    
    # Also print the formula with substituted values
    term1 = f"({d_val} - 2*{k_val})^2"
    term2_val = d_val - 2*k_val
    term2_sign = "-" if (-1)**k_val > 0 else "+"
    print(f"The calculation is: {term1} {term2_sign} {abs(term2_val)} = ({term2_val})^2 - ((-1)^{k_val})*({term2_val}) = {factor_val}")

    return factor_val

if __name__ == '__main__':
    # Example values, can be changed by the user.
    # Common case in physics is d=4.
    d_value = 4 
    k_value = 2
    
    # Calculate and print the factor
    proportionality_factor = calculate_proportionality_factor(d_value, k_value)
