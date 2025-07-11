import math

def is_prime(n):
    """Checks if an integer is a prime number."""
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

def main():
    """
    Finds and verifies solutions to x^2 - y^2 = 2023 based on given conditions.
    """
    N = 2023
    valid_pairs_count = 0
    
    print(f"Solving x^2 - y^2 = {N} for coprime integers x>0, y>0.")
    print("This requires finding factors a=x+y and b=x-y where a and b are both prime.\n")

    # Find factor pairs (a,b) where a*b=N and a>b>0
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            a = N // i
            b = i
            if a > b and b > 0:
                factor_pairs.append((a, b))

    # Process each factor pair
    for i, (a, b) in enumerate(factor_pairs):
        print(f"--- Checking Factor Pair {i+1}: a={a}, b={b} ---")

        # Compute x and y
        x = (a + b) / 2
        y = (a - b) / 2

        # x and y must be integers
        if x != int(x) or y != int(y):
            print(f"Result: x={x}, y={y} are not both integers. Invalid pair.\n")
            continue
        
        x, y = int(x), int(y)
        print(f"Computed potential solution: (x, y) = ({x}, {y})")

        # Verify all conditions
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        are_coprime = math.gcd(x, y) == 1

        print(f"Condition 1: x+y ({a}) is prime? {a_is_prime}")
        print(f"Condition 2: x-y ({b}) is prime? {b_is_prime}")
        print(f"Condition 3: x ({x}) and y ({y}) are coprime? {are_coprime}")
        
        # Check if all conditions are met
        if a_is_prime and b_is_prime and are_coprime:
            print("\nResult: All conditions are satisfied. This is a valid pair.")
            valid_pairs_count += 1
            print(f"Final Equation: {x}^2 - {y}^2 = ({x}+{y})*({x}-{y}) = {a}*{b} = {N}")
        else:
            print("\nResult: Not all conditions are met. This is not a valid pair.")
        print("-" * 25 + "\n")

    print("--- Summary ---")
    print(f"Total number of valid pairs (x,y) found: {valid_pairs_count}")

if __name__ == "__main__":
    main()
<<<0>>>