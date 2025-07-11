import math

def solve_probability():
    """
    Calculates the probability that the matrix is diagonalisable based on the
    derived formula P(diag) = 1 - (e^-1 * S3) / S2.
    """
    # The number of terms for the series approximation.
    # The factorial terms grow very rapidly, so the series converge quickly.
    # A value of 20 is more than sufficient for standard floating point precision.
    num_terms = 20

    # Calculate S2 = sum_{k=0 to inf} 1 / (k!)^2
    s2 = 0.0
    for k in range(num_terms):
        s2 += 1.0 / (math.factorial(k) ** 2)

    # Calculate S3 = sum_{k=1 to inf} 1 / (k!)^3
    s3 = 0.0
    # The sum for S3 starts from k=1
    for k in range(1, num_terms):
        s3 += 1.0 / (math.factorial(k) ** 3)

    # The value of e^-1
    exp_neg_1 = math.exp(-1)

    # Calculate the probability of not being diagonalisable
    prob_not_diag = (exp_neg_1 * s3) / s2

    # The probability of being diagonalisable is 1 minus the probability of not being so
    prob_diag = 1 - prob_not_diag

    print("The final probability is calculated using the formula: P = 1 - (e^-1 * S3) / S2")
    print("where S2 = sum_{k=0..inf} 1/(k!)^2 and S3 = sum_{k=1..inf} 1/(k!)^3.\n")
    print("The values used in the final equation are:")
    print(f"e^-1 = {exp_neg_1}")
    print(f"S2   = {s2}")
    print(f"S3   = {s3}\n")

    print("First, we calculate the probability of the matrix not being diagonalisable:")
    print(f"P(not diagonalisable) = ({exp_neg_1} * {s3}) / {s2}")
    print(f"P(not diagonalisable) = {prob_not_diag}\n")

    print("Finally, the probability of the matrix being diagonalisable is:")
    print(f"P(diagonalisable) = 1 - {prob_not_diag}")
    print(f"P(diagonalisable) = {prob_diag}")

solve_probability()