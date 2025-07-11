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

def find_solutions():
    """Finds and verifies solutions for x^2 - y^2 = 2023."""
    N = 2023
    valid_pairs_count = 0
    
    print(f"Starting analysis for the equation x^2 - y^2 = {N}")
    print("This factors to (x+y)(x-y) = 2023.")
    print("Let a = x+y and b = x-y. We are looking for factor pairs (a,b) of 2023.")
    print("The conditions are: x>0, y>0, gcd(x,y)=1, and both a and b must be prime.")
    print("-" * 50)

    # Find all factor pairs (a,b) of N where a > b
    factor_pairs = []
    for i in range(1, int(math.sqrt(N)) + 1):
        if N % i == 0:
            b = i
            a = N // i
            if a > b:
                factor_pairs.append((a, b))

    if not factor_pairs:
        print("No factor pairs found for N where a > b.")
    else:
        print(f"Found factor pairs (a, b) of {N} (with a > b): {factor_pairs}\n")

    # Check each factor pair against the conditions
    for a, b in factor_pairs:
        print(f"Checking pair a = {a}, b = {b}:")

        # Calculate potential x and y
        x = (a + b) // 2
        y = (a - b) // 2

        # Check conditions
        a_is_prime = is_prime(a)
        b_is_prime = is_prime(b)
        
        print(f"  -> Calculate x = ({a} + {b}) / 2 = {x}")
        print(f"  -> Calculate y = ({a} - {b}) / 2 = {y}")
        print(f"  -> Check 1: Is a = {a} prime? -> {'Yes' if a_is_prime else 'No'}")
        print(f"  -> Check 2: Is b = {b} prime? -> {'Yes' if b_is_prime else 'No'}")

        # If both a and b are not prime, no need to check further
        if not a_is_prime or not b_is_prime:
            print("  Result: This pair is invalid because both a and b must be prime.")
        else:
            # This part will only be reached if both a and b are prime
            are_coprime = math.gcd(x, y) == 1
            print(f"  -> Check 3: Are x = {x} and y = {y} coprime? -> {'Yes' if are_coprime else 'No'}")
            
            if are_coprime:
                print(f"  Result: Found a valid pair! (x, y) = ({x}, {y})")
                print(f"  Verification: {x}^2 - {y}^2 = {x**2} - {y**2} = {x**2 - y**2}")
                valid_pairs_count += 1
            else:
                print("  Result: This pair is invalid because x and y are not coprime.")

        print("-" * 50)

    print("\nFinal Summary:")
    print(f"The number of valid pairs (x,y) is {valid_pairs_count}.")
    
    # Return the count for the final answer block
    return valid_pairs_count

# Run the function and capture the final count
final_count = find_solutions()
print(f"<<<{final_count}>>>")