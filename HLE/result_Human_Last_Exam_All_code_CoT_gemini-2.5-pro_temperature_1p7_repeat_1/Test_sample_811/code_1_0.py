import math

def solve():
    """
    Calculates the probability that the given matrix is diagonalisable.
    """
    # Set a sufficiently large limit for numerical summation to approximate infinity
    max_k = 30
    
    # Pre-calculate p_k = P(Poisson(1) = k) = exp(-1) / k!
    p = [0.0] * (max_k + 1)
    e_neg_1 = math.exp(-1)
    for k in range(max_k + 1):
        p[k] = e_neg_1 / math.factorial(k)

    # Calculate the numerator sum S1 = sum_{k=1 to inf} p_k^3 / (1+p_k)
    s1 = 0.0
    for k in range(1, max_k + 1):
        s1 += p[k]**3 / (1 + p[k])

    # Calculate the denominator sum S2 = sum_{j=0 to inf} p_j^2 / (1+p_j)
    s2 = 0.0
    for j in range(max_k + 1):
        s2 += p[j]**2 / (1 + p[j])
        
    # The probability is 1 - S1 / S2
    probability = 1 - s1 / s2

    # Print the explanation and the numerical values as requested
    print("The probability of the matrix being diagonalisable is given by the formula:")
    print("P(diag) = 1 - S1 / S2")
    print("where S1 = sum_{k=1 to inf} p_k^3 / (1+p_k)")
    print("and S2 = sum_{j=0 to inf} p_j^2 / (1+p_j), with p_k = e^(-1)/k!.")
    print("\nNumerically computing these sums (approximated for large k):")
    print(f"S1 = {s1}")
    print(f"S2 = {s2}")
    print(f"\nThe probability is 1 - {s1} / {s2} = {probability}")

solve()