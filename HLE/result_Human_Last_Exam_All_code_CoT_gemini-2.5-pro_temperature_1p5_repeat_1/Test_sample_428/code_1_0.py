import sys

def solve_particle_annihilation():
    """
    Calculates the expected time until one particle remains based on the
    initial separations N1, M1, N2, M2.
    
    The user can provide these values as command-line arguments.
    Example: python <script_name>.py 10 20 15 25
    """
    
    # Set default values for N1, M1, N2, M2
    N1, M1, N2, M2 = 10, 20, 15, 25

    # Override with command-line arguments if provided
    if len(sys.argv) == 5:
        try:
            N1 = int(sys.argv[1])
            M1 = int(sys.argv[2])
            N2 = int(sys.argv[3])
            M2 = int(sys.argv[4])
        except ValueError:
            print("Warning: Invalid arguments. Using default values.", file=sys.stderr)
    
    # Rates for the two phases
    lambda_1 = 1.0
    lambda_2 = 2.0
    
    # --- Calculation for E[tau_1] ---
    # This is the expected time for the first collision (5 particles -> 3 particles)
    # The formula is (1/lambda_1) * sum_{i=2 to 4} (p_i - p_1) * (p_5 - p_i)
    # Let's define the terms based on N's and M's.
    # p1=0, p2=N1, p3=N1+M1, p4=N1+M1+N2, p5=N1+M1+N2+M2
    
    # Term for i=2: (p2-p1)*(p5-p2)
    p2_minus_p1 = N1
    p5_minus_p2 = M1 + N2 + M2
    term_i2 = p2_minus_p1 * p5_minus_p2
    
    # Term for i=3: (p3-p1)*(p5-p3)
    p3_minus_p1 = N1 + M1
    p5_minus_p3 = N2 + M2
    term_i3 = p3_minus_p1 * p5_minus_p3
    
    # Term for i=4: (p4-p1)*(p5-p4)
    p4_minus_p1 = N1 + M1 + N2
    p5_minus_p4 = M2
    term_i4 = p4_minus_p1 * p5_minus_p4
    
    E_tau1 = (term_i2 + term_i3 + term_i4) / lambda_1
    
    # --- Calculation for E[tau_2] ---
    # This is the expected time for the second collision (3 particles -> 1 particle)
    # The formula is (1/lambda_2) * (p_middle - p_first) * (p_last - p_middle)
    # We use the initial positions of particles p1, p3, p5.
    
    # The gaps are (p3-p1) and (p5-p3), which are already calculated.
    E_tau2 = (p3_minus_p1 * p5_minus_p3) / lambda_2

    # --- Total Expectation ---
    E_tau = E_tau1 + E_tau2

    print(f"For N1={N1}, M1={M1}, N2={N2}, M2={M2}:")
    print("-" * 30)

    # Print E[tau_1] calculation
    print("Phase 1: 5 -> 3 particles (rate = 1.0)")
    print(f"E[tau_1] = N1*(M1+N2+M2) + (N1+M1)*(N2+M2) + (N1+M1+N2)*M2")
    print(f"E[tau_1] = {N1}*({M1}+{N2}+{M2}) + ({N1}+{M1})*({N2}+{M2}) + ({N1}+{M1}+{N2})*{M2}")
    print(f"E[tau_1] = {term_i2} + {term_i3} + {term_i4} = {E_tau1}")
    print("-" * 30)

    # Print E[tau_2] calculation
    print("Phase 2: 3 -> 1 particle (rate = 2.0)")
    print(f"E[tau_2] = ( (N1+M1)*(N2+M2) ) / 2.0")
    print(f"E[tau_2] = ( ({N1}+{M1})*({N2}+{M2}) ) / 2.0")
    print(f"E[tau_2] = {term_i3} / 2.0 = {E_tau2}")
    print("-" * 30)
    
    # Print final result
    print("Total Expected Time E[tau] = E[tau_1] + E[tau_2]")
    print(f"E[tau] = ({N1}*({M1}+{N2}+{M2}) + ({N1}+{M1})*({N2}+{M2}) + ({N1}+{M1}+{N2})*{M2}) + (({N1}+{M1})*({N2}+{M2}))/2.0")
    print(f"E[tau] = {E_tau1} + {E_tau2} = {E_tau}")

if __name__ == '__main__':
    solve_particle_annihilation()