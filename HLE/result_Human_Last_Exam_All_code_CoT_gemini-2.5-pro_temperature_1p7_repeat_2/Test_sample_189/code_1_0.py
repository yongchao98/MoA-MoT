import math

def is_prime(n: int) -> bool:
    """
    Checks if a number is prime using an optimized trial division method.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    
    # Only need to check divisors up to the square root of n
    limit = int(math.sqrt(n))
    i = 5
    # Check in steps of 6 (i, i+2)
    while i <= limit:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_largest_special_prime():
    """
    Finds the largest prime p < 1,000,000 of the form p=4u+1, where u=4v+1, and p,u,v are all prime.
    """
    # Initialize final values
    p_final, u_final, v_final = (0, 0, 0)
    
    # Determine the search limit for v.
    # From p = 4(4v+1)+1 = 16v+5.
    # If p < 1,000,000, then 16v < 999,995, so v < 62,499.6875
    v_limit = 62499
    
    # Loop backwards from the limit to find the largest solution first.
    for v in range(v_limit, 1, -1):
        if is_prime(v):
            u = 4 * v + 1
            
            # The value of u must be less than the computer's limit
            if u >= 1000000:
                continue

            if is_prime(u):
                p = 4 * u + 1
                
                # The value of p must be less than the computer's limit
                if p >= 1000000:
                    continue
                    
                if is_prime(p):
                    # Since we are iterating downwards, the first solution found is the largest.
                    p_final, u_final, v_final = p, u, v
                    # We've found the largest, so we can exit the loop.
                    break

    if p_final > 0:
        print("Found the largest prime p with the specified properties.")
        print(f"The final equation values are p={p_final}, u={u_final}, v={v_final}")
        print("\nVerification of the equations:")
        print(f"{p_final} = 4 * {u_final} + 1")
        print(f"{u_final} = 4 * {v_final} + 1")
    else:
        print("No such prime p was found within the given constraints.")

# Run the function to find and print the solution.
find_largest_special_prime()