import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
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

def find_special_primes():
    """
    Finds the largest prime p < 999999 of the form p=4u+1 where u=4v+1,
    and u, v are also primes.
    """
    # The computer's integer limit is 999999.
    limit = 999999
    
    # From p = 16v + 5, we can determine the maximum possible value for v.
    # 16v + 5 < 999999 => v < 999994 / 16 => v < 62499.625
    v_limit = 62499
    
    # Variables to store the final result.
    result_p = 0
    result_u = 0
    result_v = 0
    
    # Search backwards from the v_limit to find the largest p first.
    for v_candidate in range(v_limit, 1, -1):
        # Condition 1: v must be a prime number.
        if is_prime(v_candidate):
            # Calculate u = 4v + 1
            u_candidate = 4 * v_candidate + 1
            
            # Condition 2: u must be a prime number.
            if is_prime(u_candidate):
                # Calculate p = 4u + 1
                p_candidate = 4 * u_candidate + 1
                
                # Condition 3: p must be a prime number.
                # The p < limit condition is already guaranteed by the v_limit.
                if is_prime(p_candidate):
                    # Since we are searching backwards, the first solution we find
                    # is the largest one.
                    result_p = p_candidate
                    result_u = u_candidate
                    result_v = v_candidate
                    break # Exit the loop as we've found our answer.

    # Print the final result, including the numbers in their equations.
    if result_p > 0:
        print(f"The largest prime p found is: {result_p}")
        print(f"The corresponding primes u and v are: u={result_u}, v={result_v}")
        print("\nThis means:")
        print(f"{result_p} = 4 * {result_u} + 1")
        print(f"{result_u} = 4 * {result_v} + 1")
        print(f"\nFinal Answer (p:u:v):")
        print(f"{result_p}:{result_u}:{result_v}")
    else:
        print("No such prime triplet was found within the given limit.")

# Execute the function to find and print the primes.
find_special_primes()