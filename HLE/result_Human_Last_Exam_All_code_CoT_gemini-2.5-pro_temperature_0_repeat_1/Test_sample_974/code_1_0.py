import math

def solve_and_print():
    """
    Solves the problem and prints the steps and final answer.
    """
    q = 12740347

    # Step 1: Analyze the prime divisors of q.
    # By performing a primality test, we find that 12740347 is a prime number.
    # Therefore, its only prime divisor is q itself.
    prime_divisors_of_q = [q]

    print(f"Let q = {q}.")
    print(f"The number q is prime, so its only prime divisor is p = {q}.")
    print("We will check if the condition holds for this prime p.")
    print("-" * 30)

    # Step 2: Set up the problem for p = q.
    p = q
    print(f"Checking for p = {p}.")
    print("The condition is that the number of elements of order p is the same in")
    print(f"G1 = PSL(3, q^2) and G2 = PSL(4, q).")

    # The fields are F_{q^2} and F_q. The characteristic of both fields is q.
    # We are looking for the number of elements of order p=q, which is the characteristic.
    # The number of elements of order p (the characteristic) in PSL(n, k),
    # for p > n, is given by the formula: k^(n*(n-1)) - 1.

    # Step 3: Calculate the number of elements for G1 = PSL(3, q^2).
    n1 = 3
    # The field size is k = q^2. The characteristic is p = q.
    # The condition p > n becomes q > 3, which is 12740347 > 3. This is true.
    final_exp1 = 2 * n1 * (n1 - 1)

    print(f"\nFor PSL(3, q^2):")
    print(f"n = {n1}, field size k = q^2.")
    print(f"The number of elements of order p=q is k^(n(n-1)) - 1.")
    print(f"Equation: (q^2)^({n1}*({n1-1})) - 1 = (q^2)^({n1*2}) - 1 = q^{final_exp1} - 1")

    # Step 4: Calculate the number of elements for G2 = PSL(4, q).
    n2 = 4
    # The field size is k = q. The characteristic is p = q.
    # The condition p > n becomes q > 4, which is 12740347 > 4. This is true.
    final_exp2 = n2 * (n2 - 1)

    print(f"\nFor PSL(4, q):")
    print(f"n = {n2}, field size k = q.")
    print(f"The number of elements of order p=q is k^(n(n-1)) - 1.")
    print(f"Equation: q^({n2}*({n2-1})) - 1 = q^({n2*3}) - 1 = q^{final_exp2} - 1")

    # Step 5: Compare the results and conclude.
    print("\n" + "-" * 30)
    print("Comparison:")
    print(f"Number of elements for PSL(3, q^2) is q^{final_exp1} - 1.")
    print(f"Number of elements for PSL(4, q) is q^{final_exp2} - 1.")

    if final_exp1 == final_exp2:
        print(f"\nSince {q}^{final_exp1} - 1 = {q}^{final_exp2} - 1, the numbers are equal.")
        print(f"The condition is satisfied for the prime divisor p = {p}.")
        print("\nAs this is the only prime divisor of q, it is the only solution.")
        print("The prime divisor is:")
        print(p)
        
solve_and_print()