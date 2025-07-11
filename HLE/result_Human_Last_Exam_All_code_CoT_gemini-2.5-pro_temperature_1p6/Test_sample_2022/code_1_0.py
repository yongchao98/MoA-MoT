import sys

# It is not recommended to modify the function names and parameters
# The existing code is only a template, you can modify it to suit your needs.
def solve(p):
    """
    This function takes a prime p and computes the value of F(n) mod p.
    """
    
    # Based on the analysis, the value of F(n) is S(3) mod p.
    # S(3) = 7/8. We need to compute 7 * (8^(-1)) mod p.
    
    # Python's pow(base, exp, mod) is efficient for modular exponentiation.
    # To find the modular inverse of 8 mod p, we can compute 8^(p-2) mod p.
    numerator = 7
    denominator = 8
    
    inv_denominator = pow(denominator, p - 2, p)
    
    result = (numerator * inv_denominator) % p
    
    # The problem requests to output each number in the final equation.
    # Since the logic gives the same symbolic result (7/8) for both primes,
    # the structure of the output will be the same.
    print(f"For p = {p}:")
    print(f"F(N) mod p = S(3) mod p = 7/8 mod p")
    print(f"= 7 * (8^(-1) mod {p})")
    print(f"= 7 * {inv_denominator} mod {p}")
    print(f"= {numerator * inv_denominator} mod {p}")
    print(f"= {result}")
    
    return result

def main():
    """
    Main function to solve the problem for the two given primes.
    """
    p1 = 80039
    p2 = 80077
    
    # sys.stdout.reconfigure(encoding='utf-8')
    
    result1 = solve(p1)
    print("-" * 20)
    result2 = solve(p2)
    
    print("\nFinal comma-separated answer:")
    # The final output should be in the format "<<<answer>>>"
    print(f"<<<{result1},{result2}>>>")

if __name__ == "__main__":
    main()
