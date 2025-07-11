import math

def is_prime(n):
    """
    Checks if a number n is prime using an optimized trial division method.
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
    Solves the problem by finding the required prime divisors.
    """
    q_val = 12740347

    # Step 1: Analyze the number q and its prime divisors
    print(f"Step 1: Analyze the number q = {q_val}")
    if is_prime(q_val):
        print(f"The number {q_val} is a prime number.")
        print(f"Therefore, its only prime divisor is p = {q_val}.")
        prime_divisors = [q_val]
    else:
        # Fallback for composite q, though not needed in this case
        # This part of the code will not be executed for the given q.
        print(f"The number {q_val} is not prime. This case is not expected.")
        return

    # Step 2: Explain the problem setup and the formula used
    print("\nStep 2: Set up the problem and identify the relevant formula")
    print("The problem asks for prime divisors 'p' of q for which the number of elements")
    print("of order p in PSL(3, q^2) and PSL(4, q) are equal.")
    print(f"Since p = {q_val} is the only prime divisor, p is the characteristic of the finite fields F_q and F_{q^2}.")
    print("\nThe number of elements of order p (where p is the characteristic of the field)")
    print("in the group PSL(n, k) is given by the formula:")
    print("N_p(PSL(n, k)) = k^(n * (n - 1)) - 1")

    # Step 3: Apply the formula to both groups and show the calculation
    print("\nStep 3: Apply the formula to the two groups")

    # For G1 = PSL(3, q^2)
    n1 = 3
    print(f"\n1. For G1 = PSL(3, q^2):")
    print(f"   - Here, n = {n1} and the field size is k = q^2.")
    print(f"   - The number of elements of order p=q is k^(n*(n-1)) - 1")
    print(f"     = (q^2)^({n1} * ({n1} - 1)) - 1")
    print(f"     = (q^2)^({n1 * 2}) - 1")
    print(f"     = (q^2)^6 - 1")
    print(f"     = q^12 - 1")

    # For G2 = PSL(4, q)
    n2 = 4
    print(f"\n2. For G2 = PSL(4, q):")
    print(f"   - Here, n = {n2} and the field size is k = q.")
    print(f"   - The number of elements of order p=q is k^(n*(n-1)) - 1")
    print(f"     = q^({n2} * ({n2} - 1)) - 1")
    print(f"     = q^({n2 * 3}) - 1")
    print(f"     = q^12 - 1")

    # Step 4: Conclude based on the comparison
    print("\nStep 4: Compare the results and conclude")
    print("The number of elements of order p in PSL(3, q^2) is q^12 - 1.")
    print("The number of elements of order p in PSL(4, q) is q^12 - 1.")
    print("The numbers are equal.")
    print(f"This equality holds for the prime p = {q_val}.")

    print("\nTherefore, the list of prime divisors of q satisfying the condition is:")
    for p in prime_divisors:
        print(p)

if __name__ == '__main__':
    solve()