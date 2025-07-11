import math

def solve_probability():
    """
    Calculates the probability that the given matrix is diagonalisable.
    """
    # Set a limit for the summation to approximate infinity.
    # The terms decrease very rapidly, so 30 is more than sufficient.
    limit_k = 30

    # Pre-calculate factorials for efficiency
    factorials = [1.0] * (limit_k + 1)
    for i in range(1, limit_k + 1):
        factorials[i] = factorials[i-1] * i

    # Calculate p_k = e^(-1) / k! for k from 0 to limit_k
    e_inv = math.exp(-1)
    p = [e_inv / fact for fact in factorials]

    # Calculate S2 = sum(p_k^2) and S3 = sum(p_k^3) over k from 0 to infinity
    s2 = sum(pk**2 for pk in p)
    s3 = sum(pk**3 for pk in p)

    # p0 is the probability of the Poisson variable being 0
    p0 = p[0]
    p0_cubed = p0**3
    
    # S3_star is the sum of p_k^3 from k=1 to infinity
    s3_star = s3 - p0_cubed

    # Probability of the matrix being non-diagonalisable
    prob_not_diag = s3_star / s2

    # Probability of the matrix being diagonalisable
    prob_diag = 1 - prob_not_diag

    # Print the explanation and the results
    print("The probability of the matrix being diagonalisable is calculated using the formula:")
    print("P(diag) = 1 - P(not diag)")
    print("P(not diag) = (sum_{k=1 to inf} p_k^3) / (sum_{k=0 to inf} p_k^2) = (S3 - p0^3) / S2\n")
    print("Where p_k = e^(-1)/k! is the Poisson(1) probability.")
    
    print("The required numerical values are:")
    print(f"S2 = sum p_k^2 = {s2:.15f}")
    print(f"S3 = sum p_k^3 = {s3:.15f}")
    print(f"p0^3 = (e^-1)^3 = {p0_cubed:.15f}\n")
    
    print("Plugging in the values into the formula for the probability of non-diagonalisability:")
    print(f"P(not diag) = ({s3:.15f} - {p0_cubed:.15f}) / {s2:.15f}")
    print(f"P(not diag) = {s3_star:.15f} / {s2:.15f} = {prob_not_diag:.15f}\n")
    
    print("Finally, the probability of the matrix being diagonalisable is:")
    print(f"P(diag) = 1 - {prob_not_diag:.15f} = {prob_diag:.15f}")

solve_probability()