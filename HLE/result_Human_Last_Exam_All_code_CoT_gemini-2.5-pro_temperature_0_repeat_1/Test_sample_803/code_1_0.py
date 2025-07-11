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

def is_fermat_prime(n):
    """Checks if a number is a Fermat prime."""
    if not is_prime(n):
        return False
    
    # The only known Fermat primes. For larger numbers, a more general check is needed,
    # but it's highly unlikely any other Fermat primes exist.
    known_fermat_primes = {3, 5, 17, 257, 65537}
    if n in known_fermat_primes:
        return True
    
    # General check for Fermat numbers
    if n <= 2:
        return False
    x = n - 1
    # Check if x is a power of 2
    if (x & (x - 1)) != 0:
        return False
    # Check if the exponent is a power of 2
    exponent = x.bit_length() - 1
    if (exponent & (exponent - 1)) != 0:
        return False
        
    return True

def check_if_filled(q, m):
    """
    Checks if the nonabelian group of order 2*q^m is a filled group.
    This group is assumed to be the dihedral group D_{2*q^m}.
    """
    print(f"Checking the group of order 2 * {q}^{m}:")
    
    # Validate inputs
    if not isinstance(q, int) or not isinstance(m, int) or m < 1:
        print("  Error: q must be an integer and m must be a natural number (m >= 1).")
        return
    if q % 2 == 0 or not is_prime(q):
        print(f"  Error: q must be an odd prime, but received q={q}.")
        return

    group_name = f"D_{{2 * {q}^{m}}}"
    
    is_filled_group = False
    reason = ""

    if m == 1:
        is_filled_group = True
        reason = f"Groups of the form D_{{2*q}} for an odd prime q are filled."
    elif m > 1:
        if is_fermat_prime(q):
            is_filled_group = True
            reason = f"The group is of the form D_{{2*q^m}} with m > 1, and q={q} is a Fermat prime."
        else:
            is_filled_group = False
            reason = f"The group is of the form D_{{2*q^m}} with m > 1, but q={q} is not a Fermat prime."

    if is_filled_group:
        print(f"  Result: The group {group_name} is a nonabelian filled group.")
    else:
        print(f"  Result: The group {group_name} is not a nonabelian filled group.")
    print(f"  Reason: {reason}\n")


if __name__ == '__main__':
    # Example cases
    
    # Case m=1: should be filled for any odd prime q
    check_if_filled(3, 1)
    check_if_filled(7, 1)
    
    # Case m>1: should be filled only if q is a Fermat prime
    check_if_filled(5, 2)   # 5 is a Fermat prime
    check_if_filled(11, 2)  # 11 is not a Fermat prime
    check_if_filled(3, 4)   # 3 is a Fermat prime
    check_if_filled(7, 3)   # 7 is not a Fermat prime
