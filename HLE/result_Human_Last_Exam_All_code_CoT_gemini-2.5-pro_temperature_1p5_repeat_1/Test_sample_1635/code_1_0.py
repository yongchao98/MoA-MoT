import math

def get_divisors(n):
    """
    Returns a set of all divisors of n greater than 1.
    """
    if n <= 1:
        return set()
    divs = {n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return divs

def solve_sharkovsky_problem():
    """
    Solves the problem based on the logic derived from Sharkovsky's Theorem.
    
    The problem states there is a point of order 13, but no point of order 11.
    - "Order 13 exists" implies a point of LEAST PERIOD 13 exists.
    - "No order 11" implies no point of LEAST PERIOD 11 exists.

    By Sharkovsky's theorem, the absence of a least period 11 implies the
    absence of least periods for all numbers that precede 11 in the ordering:
    3, 5, 7, 9.
    So, the set of integers k for which a point of LEAST PERIOD k does NOT exist is:
    F = {3, 5, 7, 9, 11}

    A point of "order k" exists if there's a point with a least period m,
    where m is a divisor of k (and m > 1).
    Therefore, there is NO point of order k if and only if for ALL divisors m > 1
    of k, there is NO point of least period m.
    
    This means k is in the set S if all of its divisors (>1) are in F.
    """

    # F is the set of integers k for which a point of least period k does not exist.
    F = {3, 5, 7, 9, 11}
    
    S = set()
    
    # We check numbers up to a reasonable limit to find the pattern.
    # The logic ensures that no k > 11 will be in S, except for k=1.
    limit = 100 
    
    for k in range(1, limit):
        if k == 1:
            # k=1 is special. A point of "order 1" means f(x)=x.
            # But the definition of order k requires f(x) != x.
            # So, no point of order 1 can exist by definition.
            S.add(1)
            continue
            
        divs_of_k = get_divisors(k)
        
        # Check if all divisors of k (greater than 1) are in F.
        if all(d in F for d in divs_of_k):
            S.add(k)
            
    print("Based on the problem conditions and Sharkovsky's theorem:")
    # The final equation is determining the elements of S
    print(f"S = {{k : for all m|k and m>1, m is in {sorted(list(F))}}}")
    
    sorted_S = sorted(list(S))
    
    print("\nThe numbers in set S are:")
    for number in sorted_S:
        print(number)

    print(f"\nThe full set S is: {sorted_S}")
    print(f"\nThe cardinality of S is {len(S)}.")

solve_sharkovsky_problem()
<<<6>>>