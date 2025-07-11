import math

def get_cf_coeffs(p, q):
    """Computes the continued fraction coefficients for the rational number p/q."""
    coeffs = []
    if q == 0:
        raise ValueError("Denominator cannot be zero.")
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def gcd_three(a, b, c):
    """Computes the greatest common divisor of three integers."""
    return math.gcd(math.gcd(a, b), c)

def main():
    """
    Computes the continued fraction associated with a generalized Markov number m_{p/q}
    and the quadratic equation it satisfies.
    """
    p, q = 4, 7
    print(f"The rational number is {p}/{q}.\n")

    # Step 1: Get continued fraction coefficients of the rational number
    cf_coeffs = get_cf_coeffs(p, q)
    cf_str = f"[{cf_coeffs[0]}; " + ", ".join(map(str, cf_coeffs[1:])) + "]"
    print(f"Its continued fraction is {cf_str}.\n")

    # Step 2: Extract the coefficients for construction
    # We use the part of the CF after the integer part.
    c = cf_coeffs[1:]
    print(f"The sequence of coefficients used is c = {c}.\n")

    # Step 3: Construct the palindromic period
    period = c + c[::-1]
    period_str = ", ".join(map(str, period))
    print(f"The palindromic period is constructed as c + reversed(c) = {period}.\n")

    # Step 4: Define the associated continued fraction
    print(f"The associated continued fraction is x = [{period_str}, ...].\n")
    
    # Step 5: Determine the quadratic equation for x = [period]
    # x = (p_k * x + p_{k-1}) / (q_k * x + q_{k-1})
    # This leads to q_k * x^2 + (q_{k-1} - p_k) * x - p_{k-1} = 0

    # Calculate convergents p_k, q_k, p_{k-1}, q_{k-1}
    # Initialize for convergents calculation p_n/q_n
    p_prev, p_curr = 0, 1
    q_prev, q_curr = 1, 0

    for a in period:
        p_next = a * p_curr + p_prev
        q_next = a * q_curr + q_prev
        p_prev, p_curr = p_curr, p_next
        q_prev, q_curr = q_curr, q_next
    
    p_k, p_km1 = p_curr, p_prev
    q_k, q_km1 = q_curr, q_prev

    A = q_k
    B = q_km1 - p_k
    C = -p_km1

    # Simplify the equation by dividing by GCD
    common_divisor = gcd_three(abs(A), abs(B), abs(C))
    A //= common_divisor
    B //= common_divisor
    C //= common_divisor
    
    print("This continued fraction is a root of the quadratic equation:")
    # The required output format, printing each number in the equation
    print(f"{A}x^2 + ({B})x + ({C}) = 0")

if __name__ == '__main__':
    main()