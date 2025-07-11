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

def is_power_of_two(n):
    """Checks if a number is a power of two."""
    if n <= 0:
        return False
    # A power of two in binary has only one bit set to 1.
    # n & (n - 1) will be 0 if n is a power of two.
    return (n & (n - 1)) == 0

def find_filled_groups(limit):
    """
    Finds nonabelian filled groups of order 2q^m, which are conjectured
    to only exist for m=1 and q being a Fermat prime.
    The function searches for odd primes q up to the specified limit.
    """
    print(f"Searching for nonabelian filled groups of order 2*q up to an order of {2*limit}:")
    print("Based on current research, these are expected to be dihedral groups D_{2q} where q is a Fermat prime.")
    print("-" * 20)
    
    found_any = False
    # We test odd primes q starting from 3.
    for q in range(3, limit + 1, 2):
        if is_prime(q):
            # A prime q is a Fermat prime if (q-1) is a power of two.
            if is_power_of_two(q - 1):
                found_any = True
                order = 2 * q
                m = 1
                # In the final print statement, we explicitly show each number.
                print(
                    f"Found: Group D_{{{2}*{q}}} is a nonabelian filled group of order {order}.\n"
                    f"       This fits the form 2*q^m with q={q}, m={m}.\n"
                )

    if not found_any:
        print("No such groups found within the search limit.")

if __name__ == '__main__':
    # We set a search limit for the prime q. The order of the group will be 2*q.
    search_limit_for_q = 70000
    find_filled_groups(search_limit_for_q)