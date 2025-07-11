import numpy as np

def solve():
    """
    This script demonstrates that for k=2 and d>=3, one can construct
    two valid non-isotropic covariance matrices C1 and C2 for which a
    recurrent strategy exists. This is shown by constructing them such that
    the origin is in the convex hull of their traceless parts (C'_1, C'_2).

    This serves as a counterexample for k=2, supporting the conclusion
    that the maximal k for which recurrence can never be guaranteed is 1.
    """
    # Dimension
    d = 3

    # Construct a non-zero traceless symmetric matrix E
    E = np.zeros((d, d))
    E[0, 0] = 1
    E[1, 1] = -1

    # Choose lambda > max(|eigenvalues(E)|) to ensure C1 and C2 are positive definite
    # Eigenvalues of E are 1, -1, 0. Let's pick lambda = 2.
    lam = 2.0

    # Construct the two covariance matrices
    C1 = lam * np.eye(d) + E
    C2 = lam * np.eye(d) - E

    # Check for positive definiteness
    eigenvals_C1 = np.linalg.eigvalsh(C1)
    eigenvals_C2 = np.linalg.eigvalsh(C2)

    # Calculate the traceless parts C'_i = C_i - (tr(C_i)/d) * I
    trace_C1 = np.trace(C1)
    C1_prime = C1 - (trace_C1 / d) * np.eye(d)

    trace_C2 = np.trace(C2)
    C2_prime = C2 - (trace_C2 / d) * np.eye(d)

    print(f"Demonstration for k=2, d={d}\n")
    print(f"Constructed covariance matrix C1:\n{C1}")
    print(f"Is C1 positive definite? {'Yes' if np.all(eigenvals_C1 > 0) else 'No'}")
    print(f"Is C1 isotropic (a multiple of identity)? {'No'}\n") # By construction

    print(f"Constructed covariance matrix C2:\n{C2}")
    print(f"Is C2 positive definite? {'Yes' if np.all(eigenvals_C2 > 0) else 'No'}")
    print(f"Is C2 isotropic (a multiple of identity)? {'No'}\n") # By construction

    print("For this pair of measures, a recurrent strategy exists because they are")
    print("non-degenerate and the convex hull of their traceless parts contains the origin.\n")

    print("Traceless part C'_1:")
    print(np.round(C1_prime, 5))
    print("\nTraceless part C'_2:")
    print(np.round(C2_prime, 5))
    
    sum_of_primes = C1_prime + C2_prime

    print("\nSince C'_2 = -C'_1, their sum is zero, and 0 is in their convex hull.")
    print("We demonstrate this with the equation 0.5 * C'_1 + 0.5 * C'_2 = 0")
    print("And showing that C'_1 + C'_2 equals the zero matrix:")

    # Prepare strings for the equation with each number printed
    c1p_lines = str(np.round(C1_prime, 5)).split('\n')
    c2p_lines = str(np.round(C2_prime, 5)).split('\n')
    sum_lines = str(np.round(sum_of_primes, 5)).split('\n')
    
    print("\n" + "="*20)
    print("Final Equation: C'_1 + C'_2 = 0")
    print("="*20)
    for i in range(len(c1p_lines)):
        op = '+' if i == d // 2 else ' '
        eq = '=' if i == d // 2 else ' '
        print(f"{c1p_lines[i]} {op} {c2p_lines[i]} {eq} {sum_lines[i]}")
    print("="*20)


solve()