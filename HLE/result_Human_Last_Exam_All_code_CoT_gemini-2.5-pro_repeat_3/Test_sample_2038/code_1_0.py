import snappy
import math
from collections import defaultdict

def is_prime(n):
    """
    Checks if a number n is prime.
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

def count_special_knots(max_crossings):
    """
    Counts 2-bridge knots up to max_crossings with composite determinants.
    """
    counts_by_crossing = defaultdict(int)

    # There are no knots with 1 or 2 crossings. The unknot has determinant 1, which is not composite.
    # We iterate from crossing number 3 up to max_crossings.
    for c in range(3, max_crossings + 1):
        count_for_c = 0
        # HTLinkExteriors provides the complete census of prime knots up to 15 crossings.
        # It lists one representative for a knot and its mirror image.
        for M in snappy.HTLinkExteriors(crossings=c):
            # We are interested in knots, which have 1 component.
            if M.num_components() == 1:
                # is_two_bridge() returns (p,q) for a 2-bridge knot K(p,q) or False.
                # The determinant of K(p,q) is p.
                info = M.is_two_bridge()
                if info:
                    p, q = info
                    # A number is composite if it's not prime and greater than 1.
                    # Our primality test handles n > 1 correctly.
                    # We only need to check if p is not prime, as p for a knot is always >= 3.
                    if not is_prime(p):
                        count_for_c += 1
        counts_by_crossing[c] = count_for_c
    
    return counts_by_crossing

def main():
    """
    Main function to execute the plan and print the result.
    """
    max_crossings = 13
    counts = count_special_knots(max_crossings)
    
    equation_parts = []
    # Loop from 3 to 13 to construct the equation string, matching the knot crossing numbers.
    for c in range(3, max_crossings + 1):
        equation_parts.append(str(counts.get(c, 0)))
    
    total_count = sum(counts.values())
    
    equation_str = " + ".join(equation_parts)
    
    print(f"The number of such knots for each crossing number from 3 to 13 are:")
    print(f"{equation_str} = {total_count}")

if __name__ == '__main__':
    main()
