import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_coprime_pairs_from_equation(n):
    """
    Finds all distinct pairs of coprime integers (x, y) for x^2 - y^2 = n,
    where x+y and x-y are also prime.
    """
    print(f"Finding solutions for x^2 - y^2 = {n}")
    print("This means (x+y)(x-y) = {n}.\n")
    
    factors = []
    # Find all factors of n
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            factors.append(i)
            if i*i != n:
                factors.append(n//i)
    factors.sort()
    
    # Create factor pairs (a, b) where a*b=n and a>b
    factor_pairs = []
    num_factors = len(factors)
    for i in range(num_factors // 2):
        a = factors[num_factors - 1 - i]
        b = factors[i]
        factor_pairs.append((a, b))
        
    valid_pairs_count = 0
    
    for a, b in factor_pairs:
        print(f"Considering factors a = x+y = {a}, b = x-y = {b}")
        
        # Calculate x and y
        # Since n is odd, all its factors are odd.
        # So a and b are odd, making a+b and a-b even.
        # Thus, x and y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"Resulting integer pair (x,y) = ({x}, {y})")
        print(f"Verifying the equation: {x}^2 - {y}^2 = {x**2} - {y**2} = {x**2 - y**2}")
        
        # Verify all conditions
        coprime_check = math.gcd(x, y) == 1
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"Checking conditions for (x,y) = ({x}, {y}):")
        print(f"  1. Are x and y coprime? (gcd({x}, {y}) = 1) -> {coprime_check}")
        print(f"  2. Is x+y = {a} a prime number? -> {a_is_prime}")
        print(f"  3. Is x-y = {b} a prime number? -> {b_is_prime}")
        
        if coprime_check and a_is_prime and b_is_prime:
            print(">>> All conditions met. This is a valid pair.\n")
            valid_pairs_count += 1
        else:
            print(">>> Conditions not fully met. This pair is not valid.\n")
        print("-" * 30)
        
    print(f"\nFinal Result: Found {valid_pairs_count} valid pairs.")
    return valid_pairs_count

# The number from the problem
N = 2023
num_valid_pairs = find_coprime_pairs_from_equation(N)

# The final answer in the required format
# print(f"\n<<<{num_valid_pairs}>>>")