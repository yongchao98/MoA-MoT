import sympy

def solve_particle_annihilation():
    """
    Calculates the expected time for 5 particles to annihilate down to 1.
    """
    # Define symbolic variables for the initial gap sizes
    N1, M1, N2, M2 = sympy.symbols('N1 M1 N2 M2')
    
    # Phase 1: 5 particles -> 3 particles
    # Rate lambda_1 = 1
    # Gaps are n_1=N1, n_2=M1, n_3=N2, n_4=M2
    gaps_phase1 = [N1, M1, N2, M2]
    
    # Using the formula E[tau_1] = (1/(2*lambda1)) * sum_{1<=i<=j<=k-1} n_i * n_j
    # Here k=5, lambda1=1
    sum_term_1 = 0
    for i in range(4):
        for j in range(i, 4):
            sum_term_1 += gaps_phase1[i] * gaps_phase1[j]
            
    E_tau1 = sum_term_1 / 2
    
    # Phase 2: 3 particles -> 1 particle
    # Rate lambda_2 = 2
    # We hypothesize that the relevant initial gaps for this phase are the
    # distances between particles 1, 3, and 5.
    g1 = N1 + M1
    g2 = N2 + M2
    
    # Using the formula E[tau_2] = (g1^2 + g2^2) / lambda2
    # Here lambda2 = 2
    E_tau2 = (g1**2 + g2**2) / 2
    
    # Total expected time is the sum of expectations of the two phases
    E_tau = E_tau1 + E_tau2
    
    # Expand and simplify the expression to get the final formula
    final_expr = sympy.expand(E_tau)
    
    # To conform with the output requirement "output each number in the final equation",
    # we'll construct a string representation of the formula.
    
    print("The final expression for the expectation of tau is:")
    # Using Sympy's pretty print feature
    # sympy.pprint(final_expr)
    
    # Let's extract coefficients and print them term by term.
    
    poly_expr = sympy.Poly(final_expr, N1, M1, N2, M2)
    
    terms = []
    # Loop through the monomials and their coefficients
    for monom, coeff in poly_expr.terms():
        term_str = ""
        if coeff != 1:
            # Display coefficient as a fraction if it's not an integer
            term_str += f"{sympy.pretty(sympy.nsimplify(coeff))}*"
            
        var_parts = []
        for i, p in enumerate(monom):
            if p > 0:
                var_str = [N1, M1, N2, M2][i]
                if p > 1:
                    var_parts.append(f"{var_str}**{p}")
                else:
                    var_parts.append(f"{str(var_str)}")
        term_str += "*".join(var_parts)
        terms.append(term_str)
        
    final_equation_str = "E[tau] = " + " + ".join(terms)
    
    print(final_equation_str)
    
    
    # We can also manually write out the combined and simplified expression
    # This might be more readable
    # E_tau = N1**2 + M1**2 + N2**2 + M2**2 + (3/2)*(N1*M1 + N2*M2) + (1/2)*(N1*N2 + N1*M2 + M1*N2 + M1*M2)
    # The coefficients are:
    # N1**2: 1, M1**2: 1, N2**2: 1, M2**2: 1
    # N1*M1: 3/2, N2*M2: 3/2
    # N1*N2: 1/2, N1*M2: 1/2, M1*N2: 1/2, M1*M2: 1/2
    
    
solve_particle_annihilation()
