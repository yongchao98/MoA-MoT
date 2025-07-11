import math

def main():
    """
    This program finds the limit of the sequence g_n for the given polynomial P(X)
    and demonstrates the calculation leading to the result.
    """
    
    def P(x):
        """
        Evaluates the polynomial P(X) = (X^5-1)*(X^5-X)*(X^5-X^2)*(X^5-X^3)*(X^5-X^4).
        Uses Python's arbitrary-precision integers to handle large numbers.
        """
        x_2 = x * x
        x_3 = x_2 * x
        x_4 = x_3 * x
        x_5 = x_4 * x
        
        term1 = x_5 - 1
        term2 = x_5 - x
        term3 = x_5 - x_2
        term4 = x_5 - x_3
        term5 = x_5 - x_4
        
        return term1 * term2 * term3 * term4 * term5

    def is_prime(n):
        if n <= 1: return False
        if n <= 3: return True
        if n % 2 == 0 or n % 3 == 0: return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    def find_primes(start_after, count):
        primes = []
        num = start_after + 1
        while len(primes) < count:
            if is_prime(num):
                primes.append(num)
            num += 1
        return primes

    print("To find the limit, we compute the greatest common divisor (GCD) of P(p) for several primes p.")
    print("The GCD will quickly converge to the limit.\n")

    # Let's take the first 4 primes larger than 10.
    primes = find_primes(10, 4)
    
    print(f"Step 1: Evaluate P(X) for the primes: {primes}")
    
    values = []
    for p in primes:
        val = P(p)
        values.append(val)
        print(f"For p = {p}: P({p}) = {val}")
    
    print("\nStep 2: Compute the greatest common divisor of these values.")
    
    # In Python 3.8 and below, math.gcd only takes two arguments.
    # We will compute the gcd iteratively.
    limit_candidate = values[0]
    for i in range(1, len(values)):
        limit_candidate = math.gcd(limit_candidate, values[i])

    # Show the final equation for the GCD
    gcd_all_str = f"gcd({', '.join(map(str, values))})"
    print(f"{gcd_all_str} = {limit_candidate}\n")

    print(f"The numerical calculation gives the limit as {limit_candidate}.")
    
    # Based on number theory, the limit is a product of specific prime powers.
    p1, k1 = 2, 10
    p2, k2 = 3, 2
    p3, k3 = 5, 1
    
    val1 = p1**k1
    val2 = p2**k2
    val3 = p3**k3
    
    final_limit = val1 * val2 * val3
    
    print("\nStep 3: State the theoretical result.")
    print("The limit L is the product of prime powers k_q = min_{a: gcd(a,q)=1} v_q(P(a)).")
    print(f"k_2 = 10, k_3 = 2, k_5 = 1, and k_q = 0 for all other primes q.")
    print(f"So, the limit is {p1}^{k1} * {p2}^{k2} * {p3}^{k1}:")
    print(f"{p1}^{k1} * {p2}^{k2} * {p3}^{k3} = {val1} * {val2} * {val3} = {final_limit}")

if __name__ == '__main__':
    main()