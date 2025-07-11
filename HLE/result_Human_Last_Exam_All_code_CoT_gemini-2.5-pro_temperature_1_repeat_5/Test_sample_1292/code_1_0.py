def print_field_equation():
    """
    This function prints the derived field equation for the given theory of gravity.
    The derivation is based on the variational principle applied to the action
    of Symmetric Teleparallel Gravity. The final form is taken from established
    literature and corresponds to one of the multiple-choice options.
    The output explicitly shows the numerical coefficients in the equation.
    """

    # The equation terms are constructed as strings.
    # The term 'partial_alpha(...)' denotes the partial derivative with respect to x^alpha.
    # The indices mu, nu, alpha, beta are spacetime indices.
    
    coeff1 = "-2"
    term1 = "/sqrt(-g) * partial_alpha(sqrt(-g) * P^alpha_{mu,nu})"
    
    coeff2 = "- 2"
    term2 = "* P_{mu,alpha,beta} * Q_nu^{alpha,beta}"

    # The coefficient of the third term is implicitly +1.
    coeff3 = "+ 1"
    term3 = "* Q^{alpha,beta}_mu * P_{alpha,beta,nu}"
    
    coeff4 = "- 1/2"
    term4 = "* Q * g_{mu,nu}"

    # The right-hand side of the equation.
    rhs_coeff = "8"
    rhs_term = "* pi*G/c^4 * T_{mu,nu}"

    # Print the full equation in a readable format.
    print("The derived field equation is:")
    print(f"({coeff1}) {term1} ({coeff2}) {term2} ({coeff3}) {term3} ({coeff4}) {term4} = ({rhs_coeff}) {rhs_term}")

if __name__ == "__main__":
    print_field_equation()