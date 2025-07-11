def solve():
    """
    This function calculates the value of F(n) for two prime numbers.
    The problem simplifies to calculating F(3) mod p for both primes.
    F(n) follows the recurrence F(n) = F(n-1) - (1/8) * F(n-3), with F(0)=F(1)=F(2)=1.
    Therefore, F(3) = F(2) - (1/8) * F(0) = 1 - 1/8.
    """
    
    # Prime numbers
    p1 = 80039
    p2 = 80077
    
    # Calculate modular inverse of 8 for both primes
    # using Fermat's Little Theorem: 8^(p-2) mod p
    # Or more directly with pow(8, -1, p) for Python 3.8+
    
    # For p1 = 80039
    inv8_p1 = pow(8, -1, p1)
    # Calculate F(3) for p1
    # F(3) = 1 - 8^(-1) mod p1
    f3_p1 = (1 - inv8_p1 + p1) % p1

    # For p2 = 80077
    inv8_p2 = pow(8, -1, p2)
    # Calculate F(3) for p2
    # F(3) = 1 - 8^(-1) mod p2
    f3_p2 = (1 - inv8_p2 + p2) % p2
    
    n_poly = 'p^5+2p^4-19p^3-3p^2+16p+6'

    # Print results separated by comma
    print(f"{f3_p1},{f3_p2}")

solve()
