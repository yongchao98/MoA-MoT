def print_field_equation():
    """
    Prints the derived field equation for Symmetric Teleparallel Gravity.
    """
    
    # Define the components of the equation as strings for formatting
    term1 = "-2/sqrt(-g) * d_alpha(sqrt(-g) * P^alpha_{munu})"
    term2 = "- 2*P_{mu,alpha,beta} * Q_nu^{alpha,beta}"
    term3 = "+ Q^{alpha,beta}_mu * P_{alpha,beta,nu}"
    term4 = "- 1/2 * Q * g_{munu}"
    rhs = "= (8*pi*G/c^4) * T_{munu}"

    # Print the equation term by term for clarity
    print("The derived field equation is:")
    print(f"Term 1 (Superpotential derivative): {term1}")
    print(f"Term 2 (P*Q term):                 {term2}")
    print(f"Term 3 (Q*P term):                 {term3}")
    print(f"Term 4 (Non-metricity scalar):     {term4}")
    print(f"Right-hand side (Energy-momentum): {rhs}")
    
    print("\nFinal Equation:")
    final_equation = f"{term1.replace('d_alpha', '∂_α').replace('^alpha', '^α').replace('_{munu}', '_{μν}')} " \
                     f"{term2.replace('_mu,alpha,beta', '_{μαβ}').replace('_nu^{alpha,beta}', '_ν^{αβ}')} " \
                     f"{term3.replace('^{alpha,beta}_mu', '^{αβ}_μ').replace('_{alpha,beta,nu}', '_{αβν}')} " \
                     f"{term4.replace('_{munu}', '_{μν}')} " \
                     f"{rhs.replace('_{munu}', '_{μν}')}"
    
    # A slightly more readable version using unicode
    pretty_equation = (f"-\u2202\u03b1(\u221a-g P\u03b1\u03bc\u03bd) / (\u221a-g) - 2P\u03bc\u03b1\u03b2 Q\u03bd\u03b1\u03b2 "
                       f"+ Q\u03b1\u03b2\u03bc P\u03b1\u03b2\u03bd - \u00bdQg\u03bc\u03bd = (8\u03c0G/c\u2074)T\u03bc\u03bd")
    
    # Let's print the explicit form requested by the user prompt
    print('Equation: -2 * (1/sqrt(-g)) * partial_alpha(sqrt(-g) * P^alpha_mu_nu) - 2 * P_mu_alpha_beta * Q_nu^alpha_beta + Q^alpha_beta_mu * P_alpha_beta_nu - (1/2) * Q * g_mu_nu = (8*pi*G/c^4) * T_mu_nu')
    print('Each term in the final equation is:')
    print('-g^{\rho\sigma} \, \partial_{\alpha} g_{\rho\sigma} \, P^\alpha_{\mu\nu} : Not Present')
    print(f'-2 / sqrt(-g) * partial_alpha(sqrt(-g) * P^alpha_mu_nu)')
    print(f'-2 * P_mu_alpha_beta * Q_nu^alpha_beta')
    print(f'+ Q^alpha_beta_mu * P_alpha_beta_nu')
    print(f'- (1/2) * Q * g_mu_nu')
    print(f'RHS: (8 * pi * G / c^4) * T_mu_nu')
    

print_field_equation()