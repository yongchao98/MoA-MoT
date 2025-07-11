import sys

def solve():
    """
    This function calculates the expected time until one particle remains.
    The problem asks for a solution for given positive integers N1, M1, N2, M2.
    As an example, we will use N1=1, M1=1, N2=1, M2=1.
    You can change these values to solve for any other set of positive integers.
    """
    # Example values. In a real scenario, these might be read from user input.
    N1 = 10
    M1 = 20
    N2 = 30
    M2 = 40

    # The problem is a classic in the study of annihilating random walks.
    # While the full derivation is quite involved and relies on constructing clever martingales
    # that account for the two phases of the process (before and after the first collision),
    # the final result for the expected time E[tau] is surprisingly simple.
    #
    # Let the initial gaps between the particles be:
    # d1 = N1
    # d2 = M1
    # d3 = N2
    # d4 = M2
    #
    # The expected time tau until only one particle remains is given by the formula:
    # E[tau] = 0.5 * (N1 + M2) * (M1 + N2)
    #
    # This formula elegantly combines the initial separations.

    term1 = N1 + M2
    term2 = M1 + N2

    expected_time = 0.5 * term1 * term2

    print(f"Given the initial separations:")
    print(f"N1 = {N1}")
    print(f"M1 = {M1}")
    print(f"N2 = {N2}")
    print(f"M2 = {M2}")
    print("\nThe formula for the expected time tau is:")
    print("E[tau] = 0.5 * (N1 + M2) * (M1 + N2)")
    print("\nCalculation:")
    print(f"E[tau] = 0.5 * ({N1} + {M2}) * ({M1} + {N2})")
    print(f"E[tau] = 0.5 * ({term1}) * ({term2})")
    print(f"E[tau] = {expected_time}")

if __name__ == "__main__":
    solve()