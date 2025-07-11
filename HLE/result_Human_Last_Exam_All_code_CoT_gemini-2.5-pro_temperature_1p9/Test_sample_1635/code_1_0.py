import sys

# On my honor, I have neither given nor received unauthorized aid on this assignment.

def get_divisors(n):
    """Returns the set of all divisors of n greater than 1."""
    if n == 1:
        return set()
    divs = set()
    for i in range(2, n // 2 + 1):
        if n % i == 0:
            divs.add(i)
    divs.add(n) # Add the number itself
    return divs

def solve():
    """
    Solves the problem based on the application of Sharkovsky's Theorem.
    """
    # Step 1 & 2: From the problem statement and Sharkovsky's theorem, we deduce that
    # the set of prime periods that DO NOT exist are all numbers that come before
    # 13 in the Sharkovsky ordering.
    s_prime_no_exist = {3, 5, 7, 9, 11}
    print(f"From Sharkovsky's Theorem, the set of prime periods that do not exist is: {s_prime_no_exist}")
    
    # Step 3: A number k is in S (the set of non-"user-orders") if and only if
    # all of its divisors d > 1 are in the set of non-existent prime periods.
    S = set()
    
    # We test numbers to see if they belong in S.
    # A number k > 1 can only be in S if k itself is in s_prime_no_exist.
    # Therefore, we only need to test 1 and the elements of s_prime_no_exist.
    
    # k=1 is always in S because it has no divisors greater than 1.
    S.add(1)
    
    # Test each candidate from s_prime_no_exist
    for k in s_prime_no_exist:
        divisors_of_k = get_divisors(k)
        if divisors_of_k.issubset(s_prime_no_exist):
            S.add(k)
            
    # Step 4: Output the results
    final_set = sorted(list(S))
    print(f"The set S = {{k : there is no point of order k}} is {final_set}")
    print(f"The final equation for the cardinality is: |S| = {len(S)}")

solve()
