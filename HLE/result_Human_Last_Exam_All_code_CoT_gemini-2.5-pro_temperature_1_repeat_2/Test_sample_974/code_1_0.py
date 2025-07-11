import math

def get_prime_divisors(n):
    """
    Calculates the prime divisors of a given integer n.
    """
    divisors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            divisors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        divisors.add(temp_n)
    return sorted(list(divisors))

def solve():
    """
    Solves the problem by finding prime divisors p of q for which the number
    of elements of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347

    print(f"The given number is q = {q}.")
    
    # Step 1: Find all prime divisors of q.
    prime_divisors = get_prime_divisors(q)
    print(f"The prime divisors of q are: {prime_divisors}")
    print("We need to find which of these primes p satisfy the given condition.")
    
    result_primes = []
    
    # In this case, there is only one prime divisor, p = q.
    for p in prime_divisors:
        print(f"\n--- Checking for p = {p} ---")
        
        # Step 2: Explain the theoretical background.
        # Since q is prime, its only prime divisor p is q itself.
        # The characteristic of the fields F_q and F_{q^2} is q.
        # We are looking for elements of order p=q, which is the characteristic.
        # These are called unipotent elements.
        print(f"Since p = {p} is the characteristic of the fields F_q and F_{q^2}, we need to count the number of unipotent elements.")
        print("The number of non-identity unipotent elements in PSL(n, k) is k^(n*(n-1)) - 1.")

        # Step 3: Calculate the number of elements for PSL(3, q^2)
        n1 = 3
        print(f"\nFor PSL({n1}, q^2):")
        print(f"Here, n = {n1} and the field size k is q^2.")
        print(f"The number of elements of order p is (q^2)^({n1}*({n1}-1)) - 1.")
        exp_base1 = n1 * (n1 - 1)
        final_exp1 = exp_base1 * 2
        print(f"Calculation: ({q}^2)^({exp_base1}) - 1 = {q}^{final_exp1} - 1")

        # Step 4: Calculate the number of elements for PSL(4, q)
        n2 = 4
        print(f"\nFor PSL({n2}, q):")
        print(f"Here, n = {n2} and the field size k is q.")
        print(f"The number of elements of order p is q^({n2}*({n2}-1)) - 1.")
        final_exp2 = n2 * (n2 - 1)
        print(f"Calculation: {q}^({final_exp2}) - 1")

        # Step 5: Compare the results
        print("\nComparing the two quantities:")
        print(f"PSL(3, q^2) has {q}^{final_exp1} - 1 elements of order p.")
        print(f"PSL(4, q) has {q}^{final_exp2} - 1 elements of order p.")
        
        if final_exp1 == final_exp2:
            print(f"Since {final_exp1} = {final_exp2}, the numbers are equal.")
            print(f"Therefore, p = {p} is a solution.")
            result_primes.append(p)
        else:
            print(f"Since {final_exp1} != {final_exp2}, the numbers are not equal.")

    print("\n--- Conclusion ---")
    print("The list of all prime divisors p of q for which the condition holds is:")
    print(result_primes)

if __name__ == '__main__':
    solve()