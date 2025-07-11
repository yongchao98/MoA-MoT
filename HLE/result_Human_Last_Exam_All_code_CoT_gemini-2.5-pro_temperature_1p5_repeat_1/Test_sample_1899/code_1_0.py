import decimal

def solve():
    """
    Computes the infimum and supremum of the sequence 
    E_n = |xi^n| * (a_n^2 + b_n^2 + c_n^2).
    """
    # Set precision for Decimal calculations
    decimal.getcontext().prec = 100

    # 1. Find the real root xi of f(x) = x^3 + x^2 + x - 1
    # We use Newton's method to find the root with high precision.
    # f(x) = x^3 + x^2 + x - 1
    # f'(x) = 3x^2 + 2x + 1
    xi = decimal.Decimal('0.5')  # Initial guess
    for _ in range(100):
        f_val = xi**3 + xi**2 + xi - 1
        fp_val = 3 * xi**2 + 2 * xi + 1
        if fp_val == 0:
            break
        xi = xi - f_val / fp_val

    # 2. Compute the sequence E_n for n=1 to N
    N = 200 # A sufficiently large number to observe the trend
    
    # Initial values for n=1 (P_1(x) = x)
    a, b, c = decimal.Decimal(0), decimal.Decimal(1), decimal.Decimal(0)

    E_values = []
    
    for n in range(1, N + 1):
        # Compute E_n
        norm_sq = a**2 + b**2 + c**2
        E_n = (xi**n).copy_abs() * norm_sq
        E_values.append(E_n)
        
        # Update coefficients for n+1 using the recurrence
        # a_{n+1} = c_n
        # b_{n+1} = a_n - c_n
        # c_{n+1} = b_n - c_n
        a_new = c
        b_new = a - c
        c_new = b - c
        a, b, c = a_new, b_new, c_new

    # 3. Find infimum and supremum
    infimum = min(E_values)
    supremum = max(E_values)

    print(f"Infimum: {infimum}")
    print(f"Supremum: {supremum}")

solve()
