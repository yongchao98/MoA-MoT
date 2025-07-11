import math

def is_prime(n):
    """
    Checks if a number is prime.
    """
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def get_ur(n):
    """
    Calculates the minimal order u_r(n) based on the rules established
    from the literature.
    """
    if n % 2 == 0:
        # Rule for even n
        return n // 2
    else: # n is odd
        # Rule for odd n depends on whether n is prime or composite
        if is_prime(n):
            return (n - 1) // 2
        else:
            # The only odd composite number in the range [3, 12] is 9.
            # The value for n=9 is a known special case.
            if n == 9:
                return 3
            else:
                # This case is not needed for the given problem range.
                return None

def solve():
    """
    Calculates the set of values for u_r(n) for n from 3 to 12
    and prints the result in a formatted string.
    """
    n_values = range(3, 13)
    u_r_values = [get_ur(n) for n in n_values]

    # Format the output string as a set equality
    n_str = ", ".join([f"u_r({n})" for n in n_values])
    u_r_str = ", ".join(map(str, u_r_values))
    
    print(f"{{{n_str}}} = {{{u_r_str}}}")

solve()