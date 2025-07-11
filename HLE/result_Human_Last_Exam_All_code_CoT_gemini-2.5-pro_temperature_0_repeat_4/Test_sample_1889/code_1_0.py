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

def find_factor_pairs(n):
    """Finds all factor pairs (a, b) of n such that a > b > 0."""
    pairs = []
    # Iterate from 1 up to the square root of n
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # i is a factor, so is n/i
            factor1 = i
            factor2 = n // i
            # Ensure a > b
            a = max(factor1, factor2)
            b = min(factor1, factor2)
            if a > b and b > 0:
                # Add the unique pair (a,b)
                if (a,b) not in pairs:
                    pairs.append((a, b))
    # Sort pairs by 'a' in descending order for a structured output
    pairs.sort(key=lambda p: p[0], reverse=True)
    return pairs

def solve_equation():
    """
    Finds and verifies solutions for x^2 - y^2 = 2023 based on the given conditions.
    """
    N = 2023
    print(f"The equation is x^2 - y^2 = {N}")
    print(f"This can be factored as (x+y)(x-y) = {N}")
    print("Let a = x+y and b = x-y. We need to find factors (a, b) of "
          f"{N} such that a > b > 0.")
    print("-" * 40)

    factor_pairs = find_factor_pairs(N)
    valid_pairs_count = 0
    
    if not factor_pairs:
        print(f"No factor pairs (a, b) found for {N} where a > b > 0.")
    
    for a, b in factor_pairs:
        print(f"Considering factor pair a = {a}, b = {b}")

        # Calculate x and y
        # Since N is odd, its factors a and b must be odd.
        # The sum and difference of two odd numbers are even, so x and y will be integers.
        x = (a + b) // 2
        y = (a - b) // 2
        
        print(f"  Calculating x = ({a} + {b}) / 2 = {x}")
        print(f"  Calculating y = ({a} - {b}) / 2 = {y}")

        # Verify all conditions for the pair (x, y)
        print("  Verifying conditions:")
        
        # Condition 1: x and y are coprime
        common_divisor = math.gcd(x, y)
        are_coprime = common_divisor == 1
        print(f"    1. Are x and y coprime? gcd({x}, {y}) = {common_divisor}. Result: {are_coprime}")

        # Condition 2: x+y (which is 'a') must be a prime number
        is_a_prime = is_prime(a)
        print(f"    2. Is x+y = {a} a prime number? Result: {is_a_prime}")

        # Condition 3: x-y (which is 'b') must be a prime number
        is_b_prime = is_prime(b)
        print(f"    3. Is x-y = {b} a prime number? Result: {is_b_prime}")

        if are_coprime and is_a_prime and is_b_prime:
            valid_pairs_count += 1
            print(f"  --> SUCCESS: All conditions met for (x,y) = ({x},{y}).")
        else:
            print(f"  --> FAILURE: Not all conditions were met.")
        print("-" * 40)

    print(f"\nSummary: After checking all factor pairs, the total number of valid (x, y) pairs is {valid_pairs_count}.")

if __name__ == "__main__":
    solve_equation()
    # The final answer is the count of valid pairs.
    # Based on the logic, the count will be 0.
    print("<<<0>>>")