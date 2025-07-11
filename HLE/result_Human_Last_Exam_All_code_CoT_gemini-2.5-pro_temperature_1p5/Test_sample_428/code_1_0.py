import sys

def solve():
    """
    This function calculates the expected time for the second collision.
    It reads four integer inputs N1, M1, N2, M2 from the command line.
    """
    if len(sys.argv) != 5:
        print("Usage: python your_script.py N1 M1 N2 M2")
        return

    try:
        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])
        
        if N1 <= 0 or M1 <= 0 or N2 <= 0 or M2 <= 0:
            print("Inputs N1, M1, N2, M2 must be positive integers.")
            return

    except ValueError:
        print("Invalid input. Please provide four integers.")
        return

    # Phase 1: Expected time for the first collision (5 particles -> 3 particles)
    # Rate lambda_1 = 1
    # E[tau_1] = (N1*M1 + M1*N2 + N2*M2) / (2 * lambda_1)
    e_tau1_numerator = N1 * M1 + M1 * N2 + N2 * M2
    e_tau1_denominator = 2
    e_tau1 = e_tau1_numerator / e_tau1_denominator

    # Phase 2: Expected time for the second collision (3 particles -> 1 particle)
    # Rate lambda_2 = 2
    # E[tau_rem] = E[d'_1 * d'_2] / (2 * lambda_2)
    # E[d'_1 * d'_2] = (N1+M1)*(N2+M2)
    e_tau_rem_numerator = (N1 + M1) * (N2 + M2)
    e_tau_rem_denominator = 4
    e_tau_rem = e_tau_rem_numerator / e_tau_rem_denominator

    # Total expected time E[tau]
    e_tau = e_tau1 + e_tau_rem

    print(f"N1 = {N1}")
    print(f"M1 = {M1}")
    print(f"N2 = {N2}")
    print(f"M2 = {M2}")
    print(f"E[tau] = (({N1} * {M1} + {M1} * {N2} + {N2} * {M2}) / 2) + (({N1} + {M1}) * ({N2} + {M2}) / 4)")
    print(f"E[tau] = ({e_tau1_numerator} / {e_tau1_denominator}) + ({e_tau_rem_numerator} / {e_tau_rem_denominator})")
    print(f"E[tau] = {e_tau1} + {e_tau_rem}")
    print(f"The final expectation is: {e_tau}")


if __name__ == "__main__":
    solve()
