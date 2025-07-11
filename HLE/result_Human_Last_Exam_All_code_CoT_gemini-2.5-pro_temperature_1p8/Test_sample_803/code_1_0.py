import math

def is_prime(n):
    """A simple function to check if a number is prime."""
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

def find_nonabelian_filled_groups():
    """
    Identifies and lists examples of nonabelian filled groups of order 2q^m.
    
    The identification is based on a known theorem in group theory which states that
    the only nonabelian filled groups of order 2*q^m, where q is an odd prime
    and m is a natural number, are the Dihedral groups of that order.
    
    This script provides examples for a few small values of q and m.
    """
    
    print("Based on known theorems, the nonabelian filled groups of order 2*q^m")
    print("for an odd prime q and natural number m are the Dihedral groups D_{2*q^m}.\n")
    print("Here are some examples:")
    print("=" * 40)

    # A few example values for q (odd primes) and m (natural numbers)
    odd_primes_q = [3, 5, 7]
    natural_numbers_m = [1, 2, 3]

    for q in odd_primes_q:
        # We can assert this based on our list, but a check is good practice.
        if not is_prime(q) or q == 2:
            continue
            
        for m in natural_numbers_m:
            if not isinstance(m, int) or m < 1:
                continue
            
            # Calculate the group order
            order = 2 * (q**m)
            
            # The group is the Dihedral group of this order.
            # The notation D_{order} is used for clarity.
            group_name = f"Dihedral group D_{order}"
            
            # Print the result showing each number in the equation
            print(f"For q = {q} and m = {m}:")
            print(f"  Order = 2 * {q}^{m} = {order}")
            print(f"  The filled group is the {group_name}")
            
            # Add a note for a commonly known isomorphism
            if order == 6:
                print("  (This group is also isomorphic to the Symmetric Group S_3)")
            
            print("-" * 20)

# Run the function to display the results
find_nonabelian_filled_groups()