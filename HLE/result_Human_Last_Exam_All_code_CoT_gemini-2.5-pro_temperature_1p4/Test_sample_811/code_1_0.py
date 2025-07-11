import math

def calculate_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.
    """
    # The sums converge very quickly, so 30 terms are more than sufficient
    # for standard float precision.
    N_TERMS = 30

    # Calculate S2 = sum_{j=0 to inf} 1/(j!)^2
    s2 = 0.0
    for j in range(N_TERMS):
        s2 += 1.0 / (math.factorial(j)**2)

    # Calculate S3_1 = sum_{k=1 to inf} 1/(k!)^3
    s3_from_1 = 0.0
    for k in range(1, N_TERMS):
        s3_from_1 += 1.0 / (math.factorial(k)**3)

    # Calculate the probability of not being diagonalisable
    # P(not diag) = (1/e) * S3_1 / S2
    prob_not_diag = (1.0 / math.e) * s3_from_1 / s2

    # The probability of being diagonalisable is 1 - P(not diag)
    prob_diag = 1.0 - prob_not_diag
    
    print("The probability of the matrix being diagonalisable is given by the equation:")
    print("P = 1 - (1/e) * (sum_{k=1 to inf} 1/(k!)^3) / (sum_{j=0 to inf} 1/(j!)^2)")
    print("\nLet's calculate the value of the sums and the final probability:")
    
    # Print the equation with the calculated numbers
    print(f"\nValue for S2 = sum_{j=0 to inf} 1/(j!)^2 ≈ {s2}")
    print(f"Value for S3 = sum_{k=1 to inf} 1/(k!)^3 ≈ {s3_from_1}")
    print(f"\nSubstituting these values into the equation:")
    print(f"P ≈ 1 - (1 / {math.e}) * {s3_from_1} / {s2}")
    print(f"P ≈ 1 - {prob_not_diag}")
    print(f"\nThe final probability is: {prob_diag}")

calculate_probability()