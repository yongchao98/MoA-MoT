import math

def is_prime(n):
    """
    Checks if a number is prime.
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
    Solves the problem based on the provided logic.
    """
    q = 12740347

    # Step 1: Find prime divisors of q.
    # The notation PSL(n, q) requires q to be a prime power.
    # We test if q is prime itself.
    if not is_prime(q):
        print(f"Error: q = {q} is not a prime power, the problem is ill-defined.")
        # In a real scenario, we would factorize it, but for this problem, it's prime.
        return

    print(f"The number q = {q} is a prime number.")
    # Since q is prime, its only prime divisor is p = q.
    p = q
    print(f"The only prime divisor of {q} is p = {p}.")
    print("-" * 30)

    # Step 2: Explain the logic for finding the solution.
    print("We need to check if for this prime p, the number of elements of order p")
    print(f"in PSL(3, p^2) and PSL(4, p) are equal.")
    print("\nSince p is the characteristic of the fields, elements of order p are unipotent.")
    print("The number of elements of order p in PSL(n, r) (where r is a power of p)")
    print("is given by r^(n*(n-1)) - 1, provided that p >= n.")
    print("-" * 30)

    # Step 3: Check the condition for each group.
    
    # Group 1: PSL(3, p^2)
    n1 = 3
    field_size_1_str = f"{p}^2"
    print(f"For the first group, PSL({n1}, {field_size_1_str}):")
    print(f"We check the condition p >= n, which is {p} >= {n1}.")
    if p >= n1:
        print("The condition is met.")
        num_elements_1_str = f"({field_size_1_str})^({n1}*({n1}-1)) - 1 = {p}^12 - 1"
        print(f"The number of elements of order {p} is {num_elements_1_str}.")
    else:
        print("The condition is NOT met. The formula does not apply directly.")
        # This case is not triggered for p = 12740347
    
    print("-" * 30)
    
    # Group 2: PSL(4, p)
    n2 = 4
    field_size_2_str = f"{p}"
    print(f"For the second group, PSL({n2}, {field_size_2_str}):")
    print(f"We check the condition p >= n, which is {p} >= {n2}.")
    if p >= n2:
        print("The condition is met.")
        num_elements_2_str = f"{field_size_2_str}^({n2}*({n2}-1)) - 1 = {p}^12 - 1"
        print(f"The number of elements of order {p} is {num_elements_2_str}.")
    else:
        print("The condition is NOT met. The formula does not apply directly.")
        # This case is not triggered for p = 12740347
        
    print("-" * 30)

    # Step 4: Compare and conclude.
    print("Comparing the number of elements in both groups:")
    print(f"Number in PSL(3, p^2): p^12 - 1")
    print(f"Number in PSL(4, p):   p^12 - 1")
    print("\nThe numbers are equal.")
    
    print("\nThe final equation is:")
    print(f"({p}^2)^({n1}*({n1}-1)) - 1 = {p}^({n2}*({n2}-1)) - 1")
    print(f"({p}^2)^6 - 1 = {p}^12 - 1")
    print(f"{p}^12 - 1 = {p}^12 - 1")
    
    print("\nTherefore, the condition holds for p = 12740347.")
    print("Since this is the only prime divisor of q, it is the only solution.")
    
    print("\nThe prime divisor is:")
    print(p)

solve()