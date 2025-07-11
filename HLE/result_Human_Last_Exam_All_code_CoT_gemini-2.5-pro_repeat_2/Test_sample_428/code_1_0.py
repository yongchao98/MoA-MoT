import sys

def solve():
    """
    This function calculates the expected time until only one particle remains.
    The user should provide the integer values for N1, M1, N2, M2 as command-line arguments.
    """
    if len(sys.argv) != 5:
        print("Usage: python your_script.py N1 M1 N2 M2")
        return

    try:
        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])
    except ValueError:
        print("All arguments must be integers.")
        return

    # Rate of particle movement in the first phase (5 particles)
    lambda1 = 1.0
    # Rate of particle movement in the second phase (3 particles)
    lambda2 = 2.0

    # Expected time for the first collision (from 5 to 3 particles)
    # Formula: ( (N1+M2)*(M1+N2) ) / (2*lambda1)
    E_tau1 = (N1 + M2) * (M1 + N2) / (2.0 * lambda1)

    # The expected time for the second collision (from 3 to 1 particle)
    # is a sum of two terms, derived from averaging over the 4 possible
    # collision scenarios, and includes the rate lambda2=2 in the denominator.
    # The formula from the source is:
    # 1/4 * (N1+M1)*(N2+M2) + 1/4 * (N1*M1 + N2*M2)
    # This combines the averaging and the division by lambda2=2.
    term1_tau2 = (N1 + M1) * (N2 + M2) / 4.0
    term2_tau2 = (N1 * M1 + N2 * M2) / 4.0
    E_tau2_prime = term1_tau2 + term2_tau2

    # The total expected time is the sum of the expectations of the two phases.
    total_E_tau = E_tau1 + E_tau2_prime

    print(f"Given N1={N1}, M1={M1}, N2={N2}, M2={M2}:")
    print(f"The equation for the total expected time is:")
    print(f"E[tau] = (({N1}+{M2})*({M1}+{N2}))/2 + (({N1}+{M1})*({N2}+{M2}))/4 + (({N1}*{M1})+({N2}*{M2}))/4")
    print(f"E[tau] = ({E_tau1}) + ({term1_tau2}) + ({term2_tau2})")
    print(f"The total expected time E[tau] is: {total_E_tau}")

if __name__ == "__main__":
    solve()