import sympy

def solve_annihilation_time():
    """
    Calculates the expected time for a system of 5 annihilating random walks
    to reduce to a single particle.
    """
    # Define symbolic variables for the user-provided integers
    N1 = sympy.Symbol('N_1')
    M1 = sympy.Symbol('M_1')
    N2 = sympy.Symbol('N_2')
    M2 = sympy.Symbol('M_2')

    # The initial gaps between the particles are:
    # d1 = P2 - P1 = N1
    # d2 = P3 - P2 = M1
    # d3 = P4 - P3 = N2
    # d4 = P5 - P4 = M2
    gaps = [N1, M1, N2, M2]
    gap_names = ["N_1", "M_1", "N_2", "M_2"]

    # The total "work" W to reduce 5 particles to 1 is the sum of
    # the products of all pairs of initial gaps.
    total_work = 0
    work_terms = []
    for i in range(len(gaps)):
        for j in range(i + 1, len(gaps)):
            term = gaps[i] * gaps[j]
            total_work += term
            work_terms.append(f"{gap_names[i]}*{gap_names[j]}")
    
    # This work is split equally between the two required annihilation events.
    # Work for the first collision (5 -> 3 particles)
    work_1 = total_work / 2
    # Work for the second collision (3 -> 1 particle)
    work_2 = total_work / 2

    # The rates for the two phases are given.
    rate_1 = 1
    rate_2 = 2

    # The expected time for each phase is the work divided by the rate.
    E_tau1 = work_1 / rate_1
    E_tau2 = work_2 / rate_2

    # The total expected time is the sum of the expected times for each phase.
    E_tau = E_tau1 + E_tau2
    
    # Simplify the final expression
    E_tau_simplified = sympy.simplify(E_tau)

    # Let's print out the calculation step-by-step
    print(f"The initial gaps are d_1 = N_1, d_2 = M_1, d_3 = N_2, d_4 = M_2.")
    print("The total interaction 'work' (W) is the sum of the products of all pairs of gaps:")
    print(f"W = {' + '.join(work_terms)}")
    
    print("\nThis work is divided equally between the two annihilation events.")
    print("The expected time for the first collision (E[tau_1]) is (W/2) / rate_1:")
    print(f"E[tau_1] = (({total_work}) / 2) / {rate_1} = {E_tau1}")
    
    print("\nThe expected time for the second collision (E[tau_2]) is (W/2) / rate_2:")
    print(f"E[tau_2] = (({total_work}) / 2) / {rate_2} = {E_tau2}")
    
    print("\nThe total expected time E[tau] is the sum E[tau_1] + E[tau_2]:")
    print(f"E[tau] = {E_tau1} + {E_tau2}")
    
    # Use sympy.pretty_print for a cleaner final equation
    print("\nFinal simplified expression for the total expected time E[tau]:")
    print(f"E[tau] = {E_tau_simplified}")
    

solve_annihilation_time()