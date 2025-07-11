import math

def calculate_dimension():
    """
    Calculates the largest possible dimension for R/I, which is C(10009, 10).
    """
    n = 10009
    k = 10

    # In Python 3.8+, math.comb is available and is the best way.
    # For compatibility, a manual calculation is also shown.
    try:
        result = math.comb(n, k)
    except (ImportError, AttributeError):
        if k < 0 or k > n:
            result = 0
        else:
            if k == 0 or k == n:
                result = 1
            elif k > n // 2:
                k = n - k
            
            result = 1
            for i in range(k):
                result = result * (n - i) // (i + 1)

    # Output the formula and the result as requested.
    numerator_terms = [str(i) for i in range(n, n - k, -1)]
    denominator_terms = [str(i) for i in range(k, 0, -1)]
    
    print(f"The largest possible dimension is given by the binomial coefficient C({n}, {k}).")
    print(f"C({n}, {k}) = ({' * '.join(numerator_terms)}) / ({' * '.join(denominator_terms)})")
    print(f"The final calculated value is: {result}")

calculate_dimension()