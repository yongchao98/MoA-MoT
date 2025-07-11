import sympy

def solve_particle_annihilation():
    """
    Calculates the expected time for the second collision of five particles.
    """
    # Define the symbolic variables for the initial gaps
    N1, M1, N2, M2 = sympy.symbols('N1 M1 N2 M2')

    # Initial positions of the 5 particles
    x1 = 0
    x2 = N1
    x3 = N1 + M1
    x4 = N1 + M1 + N2
    x5 = N1 + M1 + N2 + M2
    positions = [x1, x2, x3, x4, x5]

    # Step 1: Calculate E[sigma_2], the expected coalescence time for 5 particles with rate 1.
    # E[sigma_2] = (1/(2*lambda)) * Sum_{i<j} (x_j - x_i) with lambda=1.
    sum_of_distances = 0
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            sum_of_distances += (positions[j] - positions[i])

    E_sigma2 = sum_of_distances / 2
    E_sigma2_simplified = sympy.simplify(E_sigma2)

    # Step 2: Calculate E[sigma_1], the expected time for the first collision for 5 particles with rate 1.
    # The formula is E[sigma_1] = (1/lambda) * Sum_{i=1 to n-2} y_i * (Sum_{j=i+1 to n-1} y_j)
    # with n=5, lambda=1.
    # Gaps: y1=N1, y2=M1, y3=N2, y4=M2
    y = [N1, M1, N2, M2]
    E_sigma1 = 0
    for i in range(len(y) - 1):
        sum_of_later_gaps = sum(y[j] for j in range(i + 1, len(y)))
        E_sigma1 += y[i] * sum_of_later_gaps
    
    E_sigma1_simplified = sympy.simplify(E_sigma1)

    # Step 3: Combine the results using the time-change formula.
    # E[tau] = (E[sigma_1] + E[sigma_2]) / 2
    E_tau = (E_sigma1_simplified + E_sigma2_simplified) / 2
    E_tau_simplified = sympy.simplify(E_tau)

    # Print the results step-by-step
    print("The initial gaps are:")
    print(f"y1 = N1 = {N1}")
    print(f"y2 = M1 = {M1}")
    print(f"y3 = N2 = {N2}")
    print(f"y4 = M2 = {M2}\n")

    print("Step 1: Calculate E[sigma_2], the expected coalescence time in the operational time scale.")
    print(f"E[sigma_2] = (1/2) * Sum_{{1<=i<j<=5}} (x_j - x_i)")
    print(f"E[sigma_2] = {E_sigma2_simplified}\n")

    print("Step 2: Calculate E[sigma_1], the expected time of the first collision in the operational time scale.")
    print(f"E[sigma_1] = N1*(M1+N2+M2) + M1*(N2+M2) + N2*M2")
    print(f"E[sigma_1] = {E_sigma1_simplified}\n")

    print("Step 3: Calculate the final expectation E[tau].")
    print(f"E[tau] = (E[sigma_1] + E[sigma_2]) / 2")
    print(f"E[tau] = ({E_sigma1_simplified} + {E_sigma2_simplified}) / 2")
    print(f"The final expression for the expected time E[tau] is:")
    
    # To ensure the output format is a clean sum of terms
    final_expr = sympy.expand(E_tau_simplified)
    
    # Pretty print the final expression
    # This is a bit of a hack to get the desired output format without using sympy.pretty_print
    terms = []
    for term in final_expr.as_ordered_terms():
        coeff, var_part = term.as_coeff_Mul()
        if var_part == 1:
            terms.append(str(coeff))
        else:
            if coeff == 1:
                terms.append(str(var_part))
            elif coeff == -1:
                terms.append(f"-{var_part}")
            else:
                terms.append(f"{coeff}*{var_part}")
    
    final_string = " + ".join(terms).replace("+ -", "- ")
    print(final_string)

solve_particle_annihilation()