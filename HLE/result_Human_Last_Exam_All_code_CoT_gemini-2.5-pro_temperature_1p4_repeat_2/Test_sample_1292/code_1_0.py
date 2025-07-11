def print_field_equation():
    """
    This function prints the derived field equation term by term.
    The equation is of the form: J_munu = (8*pi*G/c^4) * T_munu
    We will print the terms of J_munu.
    """
    
    term1 = "-2/sqrt(-g) * d_alpha(sqrt(-g) * P^alpha_{munu})"
    term2 = "- 2 * P_{mu,alpha,beta} * Q_{nu}^{alpha,beta}"
    term3 = "+ Q^{alpha,beta}_{mu} * P_{alpha,beta,nu}"
    term4 = "- 1/2 * Q * g_{munu}"
    rhs = "= (8*pi*G/c^4) * T_{munu}"

    print("The derived field equation is:")
    print(f"({term1}) + ({term2}) + ({term3}) + ({term4}) {rhs}")
    
    print("\nSymbolically, the equation from choice B is:")
    
    # Printing each term of the final equation from choice B
    print("-2/sqrt(-g) * partial_alpha(sqrt(-g) * P^alpha_{mu,nu})", end=" ")
    print("- 2*P_{mu,alpha,beta}*Q_{nu}^{alpha,beta}", end=" ")
    print("+ Q^{alpha,beta}_{mu}*P_{alpha,beta,nu}", end=" ")
    print("- (1/2)*Q*g_{mu,nu}", end=" ")
    print("= (8*pi*G/c^4)*T_{mu,nu}")

print_field_equation()