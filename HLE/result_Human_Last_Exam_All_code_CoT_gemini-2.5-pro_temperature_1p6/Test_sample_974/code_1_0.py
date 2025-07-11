import math

def is_prime(n):
    """
    A simple primality test function to check if a number is prime.
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

def solve():
    """
    Solves the problem based on the steps outlined.
    """
    q = 12740347

    # Step 1: Explain the goal
    print(f"Let q = {q}. We need to find all prime divisors p of q for which the number of elements of order p is the same in PSL(3, q^2) and PSL(4, q).")
    
    # Step 2: Find the prime divisors of q
    print(f"\nFirst, we find the prime divisors of q = {q}.")
    if is_prime(q):
        print(f"The number q = {q} is a prime number.")
        print("Therefore, its only prime divisor is p = q itself.")
    else:
        # This part of the logic is not needed for the given q.
        print(f"The number q = {q} is not prime. Its prime divisors would need to be calculated.")
        # For the given q, this block will not execute. We proceed assuming q is prime.
        return
        
    p = q
    print(f"\nWe now check if p = {p} satisfies the condition.")
    print("The order of the elements, p = q, is the characteristic of the underlying finite fields F_q and F_{q^2}.")

    # Step 3: Apply the relevant theorem
    print("\nFor a simple group of Lie type G defined over a field of characteristic p, the number of elements of order p is given by the formula N = (|G|_p)^2 - 1.")
    print("Here, |G|_p is the order of a Sylow p-subgroup of G.")

    # Step 4: Calculation for PSL(3, q^2)
    print("\n1. For G1 = PSL(3, q^2):")
    n1 = 3
    # The Sylow q-subgroup order is |G1|_q = (q^2)^(n*(n-1)/2)
    sylow_exponent_base_power = n1 * (n1 - 1) // 2
    final_exponent_in_q = 2 * sylow_exponent_base_power
    final_number_exponent = 2 * final_exponent_in_q
    
    print(f"   The group is defined over F_(q^2), n = {n1}.")
    print(f"   The order of a Sylow q-subgroup is |G1|_q = (q^2)^({n1}({n1}-1)/2) = (q^2)^{sylow_exponent_base_power} = q^{final_exponent_in_q}.")
    print(f"   The number of elements of order q is N1 = (|G1|_q)^2 - 1 = (q^{final_exponent_in_q})^2 - 1 = q^{final_number_exponent} - 1.")

    # Step 5: Calculation for PSL(4, q)
    print("\n2. For G2 = PSL(4, q):")
    n2 = 4
    # The Sylow q-subgroup order is |G2|_q = q^(n*(n-1)/2)
    sylow_exponent_in_q = n2 * (n2 - 1) // 2
    final_number_exponent_2 = 2 * sylow_exponent_in_q
    
    print(f"   The group is defined over F_q, n = {n2}.")
    print(f"   The order of a Sylow q-subgroup is |G2|_q = q^({n2}({n2}-1)/2) = q^{sylow_exponent_in_q}.")
    print(f"   The number of elements of order q is N2 = (|G2|_q)^2 - 1 = (q^{sylow_exponent_in_q})^2 - 1 = q^{final_number_exponent_2} - 1.")
    
    # Step 6: Final Comparison and Conclusion
    print("\nComparing the two final equations:")
    print(f"Number of elements in PSL(3, q^2) = q^{final_number_exponent} - 1")
    print(f"Number of elements in PSL(4, q) = q^{final_number_exponent_2} - 1")
    
    if final_number_exponent == final_number_exponent_2:
        print("\nThe two quantities are equal.")
        print(f"Therefore, the condition holds for p = q = {q}.")
        print("\nSince this is the only prime divisor of q, it is the only solution.")
        print("\nThe list of all such prime divisors is:")
        print(p)
    else:
        print("\nThe two quantities are not equal.")

solve()