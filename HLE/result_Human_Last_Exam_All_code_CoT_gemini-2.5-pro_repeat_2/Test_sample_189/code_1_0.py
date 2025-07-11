def is_prime(n):
    """
    An efficient function to check if a number is prime.
    It handles small primes and uses the 6k +/- 1 optimization.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    # Check divisors of the form 6k +/- 1 up to sqrt(n)
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_special_prime():
    """
    Finds the largest prime p < 1,000,000 such that p = 4u + 1 and u = 4v + 1,
    with p, u, and v all being prime numbers.
    """
    # Define the upper limit for the search
    limit_p = 999999

    # Calculate the corresponding upper limits for u and v
    limit_u = (limit_p - 1) // 4
    limit_v = (limit_u - 1) // 4

    # Search backwards from the upper limit for v to find the largest result first
    for v in range(limit_v, 2, -1):
        if is_prime(v):
            # If v is prime, calculate u
            u = 4 * v + 1
            if is_prime(u):
                # If u is also prime, calculate p
                p = 4 * u + 1
                # Check if p is within the limit and is prime
                if p < limit_p and is_prime(p):
                    # Since we are searching backwards, this is the largest triplet.
                    print("Found the largest prime p, u, v triplet with p < 1,000,000.")
                    print(f"The largest prime p found is: {p}")
                    print("It is derived from primes u and v as follows:")
                    
                    # Output each number in the final equation
                    print(f"p = 4*u + 1  =>  {p} = 4*{u} + 1")
                    print(f"u = 4*v + 1  =>  {u} = 4*{v} + 1")
                    
                    print("\nAnswer in p:u:v format:")
                    print(f"{p}:{u}:{v}")
                    return

    # This part will be reached if no such triplet is found
    print("No such prime triplet was found within the specified limit.")

# Execute the function to find and print the result
find_largest_special_prime()