import sys

def solve_particle_annihilation(N1, M1, N2, M2):
    """
    Calculates the expected time until one particle remains.

    Args:
        N1 (int): A positive integer defining the initial particle separation.
        M1 (int): A positive integer defining the initial particle separation.
        N2 (int): A positive integer defining the initial particle separation.
        M2 (int): A positive integer defining the initial particle separation.
    """
    if not all(isinstance(i, int) and i > 0 for i in [N1, M1, N2, M2]):
        print("Error: N1, M1, N2, M2 must be positive integers.", file=sys.stderr)
        return

    # E[tau1] is the sum of products of non-adjacent initial gaps.
    E_tau1 = N1 * N2 + N1 * M2 + M1 * M2

    # E[tau2] is half the sum of products of adjacent initial gaps.
    E_tau2 = 0.5 * (N1 * M1 + M1 * N2 + N2 * M2)

    # The total expected time is the sum of the two expectations.
    E_tau = E_tau1 + E_tau2

    print("The problem is solved by calculating the expected durations of two phases:")
    print("Phase 1 (5 to 3 particles, rate 1): E[tau1]")
    print("Phase 2 (3 to 1 particle, rate 2): E[tau2]")
    print("\nBased on the initial gaps N1, M1, N2, M2:")
    print(f"N1 = {N1}, M1 = {M1}, N2 = {N2}, M2 = {M2}")
    
    print("\nThe expectation for the first phase is the sum of products of non-adjacent gaps:")
    print(f"E[tau1] = N1*N2 + N1*M2 + M1*M2 = {N1}*{N2} + {N1}*{M2} + {M1}*{M2} = {E_tau1}")
    
    print("\nThe expectation for the second phase is half the sum of products of adjacent gaps:")
    print(f"E[tau2] = 0.5 * (N1*M1 + M1*N2 + N2*M2) = 0.5 * ({N1}*{M1} + {M1}*{N2} + {N2}*{M2}) = {E_tau2}")
    
    print("\nThe total expected time E[tau] is the sum E[tau1] + E[tau2]:")
    print(f"E[tau] = {E_tau1} + {E_tau2} = {E_tau}")

# Example values for N1, M1, N2, M2.
# You can change these values to see the result for different initial separations.
N1 = 10
M1 = 20
N2 = 10
M2 = 30

solve_particle_annihilation(N1, M1, N2, M2)