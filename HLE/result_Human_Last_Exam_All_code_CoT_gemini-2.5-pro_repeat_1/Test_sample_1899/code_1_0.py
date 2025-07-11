import decimal

def solve_polynomial_sequence():
    """
    Calculates the infimum and supremum of the sequence 
    E_n = |xi|^n * (a_n^2 + b_n^2 + c_n^2)
    numerically.
    """
    # Set precision for decimal arithmetic
    # The coefficients grow large, and xi^n becomes small, so high precision is needed.
    prec = 100
    decimal.getcontext().prec = prec

    # 1. Find the real root xi of x^3 + x^2 + x - 1 = 0
    # We use Newton's method for high precision.
    f = lambda x: x**3 + x**2 + x - 1
    f_prime = lambda x: 3*x**2 + 2*x + 1
    
    # Initial guess
    x = decimal.Decimal('0.5')
    
    # Iterate Newton's method
    for _ in range(20):
        x = x - f(x) / f_prime(x)
    
    xi = x
    
    # 2. Initialize coefficients for n=1
    # For n=1, xi^1 = a1 + b1*xi + c1*xi^2 => (a1,b1,c1) = (0,1,0)
    a, b, c = 0, 1, 0

    # Initialize min and max E_n with the value for n=1
    xi_n = xi
    E_n = xi_n * (a**2 + b**2 + c**2)
    min_E = E_n
    max_E = E_n

    # 3. Iterate to find inf and sup
    num_iterations = 2000
    for n in range(2, num_iterations + 1):
        # Update coefficients using the recurrence relation:
        # a_new = c
        # b_new = a - c
        # c_new = b - c
        a_new = c
        b_new = a - c
        c_new = b - c
        a, b, c = a_new, b_new, c_new
        
        # Calculate xi^n
        xi_n *= xi
        
        # Calculate E_n
        E_n = xi_n * (a**2 + b**2 + c**2)
        
        # Update min and max
        if E_n < min_E:
            min_E = E_n
        if E_n > max_E:
            max_E = E_n
    
    # Set precision for final output
    decimal.getcontext().prec = 15
    print(f"The real root xi is approximately: {+xi}")
    print(f"After {num_iterations} iterations:")
    print(f"Infimum found: {+min_E}")
    print(f"Supremum found: {+max_E}")

solve_polynomial_sequence()