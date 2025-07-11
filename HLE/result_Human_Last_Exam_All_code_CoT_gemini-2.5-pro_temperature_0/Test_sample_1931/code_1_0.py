import random

def is_prime(n):
    """
    Performs the Miller-Rabin primality test, which is efficient for large numbers.
    It's a probabilistic test, but highly reliable for this purpose.
    """
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    if n < 25:
        return True
    
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Run 5 rounds of the test for high confidence.
    for _ in range(5):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def solve_riddle():
    """
    Finds the 6th prime in the sequence of concatenated digits of 'e'
    and verifies its last 6 digits.
    """
    # The digits of the irrational number 'e', starting after the decimal point.
    # We need enough digits to find the 6th prime, which is 353 digits long.
    e_digits = "27182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841675092447614606680822648001684774118537423454424371075390777449920695506334540650341217317318477629296333894947245535947762495061800782395239880232344476298412553222220550362242321463135714636811135950163483189174424333861013414423823758521023"
    
    target_prime_index = 6
    target_suffix = "521023"
    prime_count = 0
    
    print("The irrational number is 'e'.")
    print("Searching for the first 6 primes from its concatenated digits...")
    print("-" * 60)

    for i in range(1, len(e_digits) + 1):
        num_str = e_digits[:i]
        
        # Optimization: Primes (except 2, 5) can't end in an even digit or 5.
        if i > 1 and num_str[-1] in '024568':
            continue
            
        num = int(num_str)
        if is_prime(num):
            prime_count += 1
            print(f"Found prime #{prime_count} (length {len(num_str)}): {num_str[:30]}...")
            
            if prime_count == target_prime_index:
                print("-" * 60)
                print(f"Found the {target_prime_index}th prime number.")
                
                if num_str.endswith(target_suffix):
                    print(f"Its last 6 digits are indeed '{target_suffix}'.")
                    print("\nThe riddle is solved! The irrational number is 'e'.")
                    print("\nThe full prime number (the 'equation' number) is:")
                    print(num_str)
                else:
                    print("Error: The last 6 digits did not match the riddle's criteria.")
                break

if __name__ == "__main__":
    solve_riddle()