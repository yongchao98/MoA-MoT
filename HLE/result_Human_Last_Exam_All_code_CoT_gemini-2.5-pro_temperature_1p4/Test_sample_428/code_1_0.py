def solve_particle_annihilation():
    """
    Calculates the expected time for the second collision of 5 annihilating particles.
    """
    # Positive integer inputs for the problem.
    # You can change these values to solve for different initial configurations.
    N1 = 1
    M1 = 1
    N2 = 1
    M2 = 1

    # Stage 1: Expected time from 5 to 3 particles, with rate lambda_1 = 1.
    
    # The sum of products of gaps for all 10 triplets of particles
    # simplifies to the following expression:
    # S = 3*N1*M1 + 4*N1*N2 + 3*N1*M2 + 4*M1*N2 + 3*M1*M2 + 3*N2*M2
    # The expected time for the first collision is S / (1 * (5-1)*(5-2)) = S / 6
    
    term_N1M1 = 3 * N1 * M1
    term_N1N2 = 4 * N1 * N2
    term_N1M2 = 3 * N1 * M2
    term_M1N2 = 4 * M1 * N2
    term_M1M2 = 3 * M1 * M2
    term_N2M2 = 3 * N2 * M2
    
    sum_of_prods = term_N1M1 + term_N1N2 + term_N1M2 + term_M1N2 + term_M1M2 + term_N2M2
    
    lambda_1 = 1
    e_tau1 = sum_of_prods / (lambda_1 * (5 - 1) * (5 - 2))
    
    # Stage 2: Expected time from 3 to 1 particle, with rate lambda_2 = 2.
    lambda_2 = 2
    d = [N1, M1, N2, M2]
    
    # Probabilities of which gap closes first (p_i is proportional to 1/d_i)
    sum_inv_d = sum(1.0/val for val in d)
    probs = [(1.0/val) / sum_inv_d for val in d]
    
    # Expected duration of stage 2, conditioned on which pair collides first.
    # We use the initial gaps of the surviving particles as a proxy.
    exp_duration_cond = [
        (N2 * M2) / lambda_2,
        ((N1 + M1 + N2) * M2) / lambda_2,
        (N1 * (M1 + N2 + M2)) / lambda_2,
        (N1 * M1) / lambda_2
    ]
    
    e_tau2 = sum(p * E for p, E in zip(probs, exp_duration_cond))
    
    e_tau = e_tau1 + e_tau2

    # Outputting the step-by-step calculation with the given numbers
    print(f"Given N1={N1}, M1={M1}, N2={N2}, M2={M2}:")
    print("\n--- Step 1: Calculate E[tau_1] (time for first collision, rate=1) ---")
    print(f"E[tau_1] = (1/6) * (3*N1*M1 + 4*N1*N2 + 3*N1*M2 + 4*M1*N2 + 3*M1*M2 + 3*N2*M2)")
    print(f"E[tau_1] = (1/6) * (3*{N1}*{M1} + 4*{N1}*{N2} + 3*{N1}*{M2} + 4*{M1}*{N2} + 3*{M1}*{M2} + 3*{N2}*{M2})")
    print(f"E[tau_1] = (1/6) * ({term_N1M1} + {term_N1N2} + {term_N1M2} + {term_M1N2} + {term_M1M2} + {term_N2M2})")
    print(f"E[tau_1] = (1/6) * {sum_of_prods} = {e_tau1}")

    print("\n--- Step 2: Calculate E[tau_2] (time for second collision, rate=2) ---")
    print("This requires averaging over the 4 possibilities for the first collision.")
    print(f"Probabilities (p1, p2, p3, p4) = ({probs[0]:.4f}, {probs[1]:.4f}, {probs[2]:.4f}, {probs[3]:.4f})")
    print("Conditional expected durations E_i (rate=2):")
    print(f"E_1 = (N2*M2)/2 = ({N2}*{M2})/2 = {exp_duration_cond[0]}")
    print(f"E_2 = ((N1+M1+N2)*M2)/2 = (({N1}+{M1}+{N2})*{M2})/2 = {exp_duration_cond[1]}")
    print(f"E_3 = (N1*(M1+N2+M2))/2 = ({N1}*({M1}+{N2}+{M2}))/2 = {exp_duration_cond[2]}")
    print(f"E_4 = (N1*M1)/2 = ({N1}*{M1})/2 = {exp_duration_cond[3]}")
    print(f"E[tau_2] = p1*E_1 + p2*E_2 + p3*E_3 + p4*E_4")
    print(f"E[tau_2] = {probs[0]:.4f}*{exp_duration_cond[0]} + {probs[1]:.4f}*{exp_duration_cond[1]} + {probs[2]:.4f}*{exp_duration_cond[2]} + {probs[3]:.4f}*{exp_duration_cond[3]} = {e_tau2:.4f}")

    print("\n--- Step 3: Total Expected Time E[tau] ---")
    print(f"E[tau] = E[tau_1] + E[tau_2]")
    print(f"E[tau] = {e_tau1} + {e_tau2:.4f} = {e_tau:.4f}")

solve_particle_annihilation()