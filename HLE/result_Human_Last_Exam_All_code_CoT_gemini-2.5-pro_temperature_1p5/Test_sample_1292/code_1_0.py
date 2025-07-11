def derive_field_equation():
    """
    This function prints the derived field equation for Symmetric Teleparallel Gravity.
    The derivation is based on the variational principle applied to the given action.
    """
    
    # Define the components of the equation as strings
    term1 = "-2/sqrt(-g) * d_alpha(sqrt(-g) * P^alpha_{munu})"
    term2 = "- 2*P_{mu,alpha,beta} * Q_nu^{alpha,beta}"
    term3 = "+ Q^{alpha,beta}_mu * P_{alpha,beta,nu}"
    term4 = "- 1/2 * Q * g_{munu}"
    
    # RHS of the equation
    rhs = "= 8*pi*G/c^4 * T_{munu}"
    
    # Combine terms to form the full equation
    field_equation = f"{term1} {term2} {term3} {term4} {rhs}"

    print("The derived field equation is:")
    
    # Print each part of the equation symbolically
    print(
        f"Term 1 (derivative of superpotential): "
        f"-\u2202\u03b1(\u221A-g P\u03b1\u03bc\u03bd) * (2/\u221A-g)  "
        f"which corresponds to the code output part: {term1}\n"

        f"Term 2 (interaction term 1): "
        f"-2P\u03bc\u03b1\u03b2 Q\u03bd\u03b1\u03b2  "
        f"which corresponds to the code output part: {term2}\n"

        f"Term 3 (interaction term 2): "
        f"+Q\u03b1\u03b2\u03bc P\u03b1\u03b2\u03bd  "
        f"which corresponds to the code output part: {term3}\n"

        f"Term 4 (non-metricity scalar term): "
        f"-(1/2)Q g\u03bc\u03bd  "
        f"which corresponds to the code output part: {term4}\n"
    )

    print("Final Equation:")
    # Print the equation in a more readable format using unicode characters for symbols
    # ∂α(√-g Pαμν) * (2/√-g) - 2Pμαβ Qναβ + Qαβμ Pαβν - (1/2)Q gμν = (8πG/c⁴) Tμν
    final_equation_pretty = (
        f"-\u2202\u03b1(\u221A-g P^\u03b1\u03bc\u03bd) * (2/\u221A-g) - "
        f"2P\u03bc\u03b1\u03b2 Q^\u03b1\u03b2\u03bd + "
        f"Q^\u03b1\u03b2\u03bc P\u03b1\u03b2\u03bd - "
        f"(1/2)Q g\u03bc\u03bd = (8\u03c0G/c\u2074) T\u03bc\u03bd"
    )

    # Simplified form matching option B.
    final_equation_b = (
        f"-\frac{{2}}{{\sqrt{{-g}}}}\partial_{{\alpha}}(\sqrt{{-g}}P^\alpha_{{\mu\nu}}) - "
        f"2P_{{\mu\alpha\beta}} Q_\nu^{{\alpha\beta}} + "
        f"Q^{{\alpha\beta}}_\mu P_{{\alpha\beta\nu}} - "
        f"\\frac{{1}}{{2}}Qg_{{\mu\nu}}= "
        f"\\frac{{8\pi G}}{{c^4}} T_{{\mu\nu}}"
    )
    
    # We output each part of the final equation to show the numbers
    print(f"-2 * (1/sqrt(-g)) * d_alpha(sqrt(-g) * P^alpha_munu) - 2 * P_mu_alpha_beta * Q_nu^alpha_beta + 1 * Q^alpha_beta_mu * P_alpha_beta_nu - 1/2 * Q * g_munu = (8*pi*G/c^4) * T_munu")


derive_field_equation()