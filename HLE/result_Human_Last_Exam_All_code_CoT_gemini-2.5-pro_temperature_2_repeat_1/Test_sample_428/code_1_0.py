def solve():
    """
    Calculates the expected time tau for the second collision to happen.

    The user is expected to provide the positive integers N1, M1, N2, M2.
    These define the initial locations of the five particles as:
    0, N1, N1+M1, N1+M1+N2, N1+M1+N2+M2.

    The expected time is calculated in two parts:
    1. E[tau_1]: The expected time for the first collision (5 particles -> 3 particles).
       During this phase, all particles move with rate 1.
       E[tau_1] = (N1*M1 + M1*N2 + N2*M2) / 1.

    2. E[tau_2]: The expected time for the second collision (3 particles -> 1 particle).
       During this phase, the remaining three particles move with rate 2.
       The initial gaps for this phase depend on which particles collided first.
       Averaging over all possibilities gives a beautiful result in terms of the
       initial non-adjacent gaps.
       E[tau_2] = (N1*N2 + N1*M2 + M1*M2) / 2.

    The total expected time E[tau] is the sum of these two expectations.
    """
    try:
        n1_str = input("Enter the value of N1: ")
        m1_str = input("Enter the value of M1: ")
        n2_str = input("Enter the value of N2: ")
        m2_str = input("Enter the value of M2: ")
        N1 = int(n1_str)
        M1 = int(m1_str)
        N2 = int(n2_str)
        M2 = int(m2_str)

        if not all(x > 0 for x in [N1, M1, N2, M2]):
            print("All values N1, M1, N2, M2 must be positive integers.")
            return

        # Calculate E[tau_1] (rate = 1)
        # Sum of products of adjacent initial gaps
        e_tau1 = float(N1 * M1 + M1 * N2 + N2 * M2)

        # Calculate E[tau_2] (rate = 2)
        # Sum of products of non-adjacent initial gaps, divided by the new rate
        e_tau2 = float(N1 * N2 + N1 * M2 + M1 * M2) / 2.0

        # Total expected time
        total_expected_time = e_tau1 + e_tau2

        # Print the final equation with all numbers
        print("The expected time E[tau] is calculated as E[tau_1] + E[tau_2]")
        print("where E[tau_1] corresponds to the first collision (rate=1) and E[tau_2] to the second (rate=2).")
        print("\nBased on the initial gaps between particles (N1, M1, N2, M2):")
        
        # Breakdown of the calculation for clarity
        # Products of adjacent gaps for E[tau_1]
        term_n1_m1 = N1 * M1
        term_m1_n2 = M1 * N2
        term_n2_m2 = N2 * M2
        
        # Products of non-adjacent gaps for E[tau_2]
        term_n1_n2 = N1 * N2
        term_n1_m2 = N1 * M2
        term_m1_m2 = M1 * M2
        
        # E[tau_1] and E[tau_2] values
        val_e_tau1 = float(term_n1_m1 + term_m1_n2 + term_n2_m2)
        val_e_tau2 = float(term_n1_n2 + term_n1_m2 + term_m1_m2) / 2.0
        
        print(f"\nE[tau_1] = (N1*M1 + M1*N2 + N2*M2) / 1")
        print(f"E[tau_1] = ({N1}*{M1} + {M1}*{N2} + {N2}*{M2}) / 1")
        print(f"E[tau_1] = ({term_n1_m1} + {term_m1_n2} + {term_n2_m2}) / 1 = {val_e_tau1}")

        print(f"\nE[tau_2] = (N1*N2 + N1*M2 + M1*M2) / 2")
        print(f"E[tau_2] = ({N1}*{N2} + {N1}*{M2} + {M1}*{M2}) / 2")
        print(f"E[tau_2] = ({term_n1_n2} + {term_n1_m2} + {term_m1_m2}) / 2 = {val_e_tau2}")

        print(f"\nE[tau] = E[tau_1] + E[tau_2] = {val_e_tau1} + {val_e_tau2} = {total_expected_time}")
        print(f"\nFinal Answer: {total_expected_time}")

    except ValueError:
        print("Invalid input. Please enter positive integer values.")
    except Exception as e:
        print(f"An error occurred: {e}")

solve()
<<<E[\tau] = N_1 M_1 + M_1 N_2 + N_2 M_2 + \frac{1}{2}(N_1 N_2 + N_1 M_2 + M_1 M_2)>>>