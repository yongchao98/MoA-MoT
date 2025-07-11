import math

def get_prime_factors(n):
    """
    Returns a list of unique prime factors of a given integer n.
    """
    factors = set()
    temp = n
    # Check for factor 2
    if temp % 2 == 0:
        factors.add(2)
        while temp % 2 == 0:
            temp //= 2
    # Check for odd factors
    d = 3
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 2
    if temp > 1:
        factors.add(temp)
    return sorted(list(factors))

def find_solution_primes():
    """
    Finds and prints the list of prime divisors p of q for which the number of
    elements of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347
    
    print(f"Given q = {q}. We seek its prime divisors p that satisfy the condition.")
    
    prime_divisors = get_prime_factors(q)
    
    print(f"The prime divisors of q are: {prime_divisors}\n")

    result_primes = []
    
    # The condition for the number of elements of order p to be equal
    # simplifies to requiring p >= 5 for the relevant formulas to hold.
    for p in prime_divisors:
        print(f"--- Checking prime divisor p = {p} ---")
        
        condition_met = (p >= 5)
        
        if condition_met:
            print(f"The condition p >= 5 is met.")
            
            # For PSL(3, q^2), n=3. The number of elements of order p is (q^2)^(3*(3-1)) - 1
            n1 = 3
            exp1_base = 2
            exp1_power = n1 * (n1 - 1)
            final_exp1 = exp1_base * exp1_power
            
            # For PSL(4, q), n=4. The number of elements of order p is q^(4*(4-1)) - 1
            n2 = 4
            final_exp2 = n2 * (n2 - 1)

            print("Equation for number of elements of order p:")
            print(f"In PSL({n1}, q^{exp1_base}): (q^{exp1_base})^({n1}*({n1}-1)) - 1 = q^{final_exp1} - 1")
            print(f"In PSL({n2}, q): q^({n2}*({n2}-1)) - 1 = q^{final_exp2} - 1")
            
            if final_exp1 == final_exp2:
                print("\nThe exponents are equal, so the number of elements are equal.")
                result_primes.append(p)
            # This 'else' will not be reached as 12 == 12.
            else:
                 print("\nThe exponents are not equal.")
            
        else:
            print(f"The condition p >= 5 is NOT met.")
        print("-" * 35)

    print("\nFinal list of all prime divisors p for which the number of elements are equal:")
    for p in result_primes:
        print(p)

find_solution_primes()