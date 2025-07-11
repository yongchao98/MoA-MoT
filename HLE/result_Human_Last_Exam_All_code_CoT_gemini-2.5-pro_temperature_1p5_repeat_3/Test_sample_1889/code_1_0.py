import math

def is_prime(n):
    """
    Checks if a number n is prime.
    A number is prime if it is greater than 1 and has no divisors other than 1 and itself.
    """
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

def solve_equation():
    """
    Finds solutions to x^2 - y^2 = 2023 with the given constraints.
    """
    N = 2023
    valid_solutions_count = 0

    print(f"Solving x^2 - y^2 = {N} for positive, coprime integers (x, y).")
    print("This requires that both x+y and x-y are prime numbers.\n")
    print(f"Step 1: Factor the equation to (x+y)(x-y) = {N}.")
    print("Let a = x+y and b = x-y. We must find factor pairs (a, b) of 2023.\n")

    # Find factor pairs (a, b) of N where a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            a = N // i
            b = i
            if a > b:
                factor_pairs.append((a, b))

    print(f"Step 2: The factor pairs (a, b) of {N} where a > b are: {factor_pairs}\n")
    print("Step 3: Check if for any pair, both a and b are prime.")

    for a, b in factor_pairs:
        print(f"\n----- Checking Pair (a, b) = ({a}, {b}) -----")
        
        # Check if both a and b are prime
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"  Condition: Is a = {a} (which is x+y) a prime number? {a_is_prime}")
        print(f"  Condition: Is b = {b} (which is x-y) a prime number? {b_is_prime}")

        if a_is_prime and b_is_prime:
            # If both are prime, calculate x and y and check coprimality
            print("  Result: Both a and b are prime. Proceeding to calculate x and y.")
            x = (a + b) // 2
            y = (a - b) // 2
            
            # For distinct odd primes a,b, x and y are guaranteed to be coprime.
            # gcd((a+b)/2, (a-b)/2) = gcd(a,b) = 1.
            are_coprime = math.gcd(x, y) == 1

            print(f"  Calculated x = ({a} + {b}) / 2 = {x}")
            print(f"  Calculated y = ({a} - {b}) / 2 = {y}")
            print(f"  Condition: Are x and y coprime? gcd({x}, {y}) = {math.gcd(x, y)}. {'Yes' if are_coprime else 'No'}")
            
            if are_coprime:
                print(f"\n  SUCCESS: A valid pair is found: (x, y) = ({x}, {y})")
                print(f"  Final Equation: {x}^2 - {y}^2 = {x*x} - {y*y} = {N}")
                valid_solutions_count += 1
        else:
            print("  Result: The pair is invalid because both a and b must be prime.")

    print("\n----- Summary -----")
    if valid_solutions_count == 0:
        print("No factor pair (a, b) consists of two prime numbers.")
    
    print(f"The total number of valid pairs (x, y) is {valid_solutions_count}.")

# Run the solver
solve_equation()
<<<0>>>