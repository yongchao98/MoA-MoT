from fractions import Fraction

def calculate_expected_time(N1, M1, N2, M2):
    """
    Calculates the expected time tau for the second collision to happen.

    Args:
        N1, M1, N2, M2: Positive integers representing the initial separations.
    """
    # E[tau] = 1/3 * (N1*M1 + M1*N2 + N2*M2) + 5/6 * (N1*N2 + N1*M2 + M1*M2)

    # Calculate the coefficients as Fraction objects for precision
    c1 = Fraction(1, 3)
    c2 = Fraction(5, 6)

    # Calculate each term of the final equation
    term1 = c1 * N1 * M1
    term2 = c1 * M1 * N2
    term3 = c1 * N2 * M2
    term4 = c2 * N1 * N2
    term5 = c2 * N1 * M2
    term6 = c2 * M1 * M2

    print("The final equation for the expected time E[tau] is a sum of the following terms:")
    print(f"Term 1 (from N1*M1): {N1*M1}/3 = {float(term1)}")
    print(f"Term 2 (from M1*N2): {M1*N2}/3 = {float(term2)}")
    print(f"Term 3 (from N2*M2): {N2*M2}/3 = {float(term3)}")
    print(f"Term 4 (from N1*N2): 5*{N1*N2}/6 = {float(term4)}")
    print(f"Term 5 (from N1*M2): 5*{N1*M2}/6 = {float(term5)}")
    print(f"Term 6 (from M1*M2): 5*{M1*M2}/6 = {float(term6)}")

    # Sum all terms to get the final result
    total_expectation = term1 + term2 + term3 + term4 + term5 + term6

    print("\nTotal expected time E[tau]:")
    print(f"E[tau] = {float(total_expectation)}")

if __name__ == '__main__':
    # Example values for the positive integers N1, M1, N2, M2.
    # You can change these to any positive integer values.
    N1 = 10
    M1 = 20
    N2 = 15
    M2 = 25

    print(f"Calculating for N1={N1}, M1={M1}, N2={N2}, M2={M2}\n")
    calculate_expected_time(N1, M1, N2, M2)