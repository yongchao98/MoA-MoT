import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    It counts the positive integers up to a given integer n
    that are relatively prime to n.
    """
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve():
    """
    Finds the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    k = 6
    primes = [7, 13, 19, 31] # Prime factorization of 53599

    print(f"We want to find the number of primitive Dirichlet characters of conductor d = {d} and order k = {k}.")
    print(f"The conductor d = {d} is square-free, with prime factors: {primes}.")
    print("A character chi is primitive with conductor d if and only if it is a product of non-principal characters")
    print("mod each prime factor p of d. The order of chi is the lcm of the orders of its components.")
    
    # We need lcm(ord(chi_7), ord(chi_13), ord(chi_19), ord(chi_31)) = 6.
    # This means the orders of the components must be divisors of 6, i.e., {2, 3, 6} (since they must be > 1).
    # Also, for an order m to be possible for a character mod p, m must divide p-1.
    print("\nFor each prime p, the order of a component character mod p must divide p-1:")
    for p in primes:
        print(f"p = {p}, p-1 = {p-1}. Since {k} divides {p-1}, characters of order {k} exist.")

    # Using the Principle of Inclusion-Exclusion.
    # Total = (char with order dividing 6) - (char with order dividing 3) - (char with order dividing 2)
    
    # 1. Count combinations where component orders divide 6
    # For each prime, possible primitive orders are 2, 3, 6.
    # Number of choices = phi(2) + phi(3) + phi(6) = 1 + 2 + 2 = 5
    num_choices_div_6 = phi(2) + phi(3) + phi(6)
    total_div_6 = num_choices_div_6 ** len(primes)
    print(f"\nStep 1: Count combinations where component orders are in {{2, 3, 6}}.")
    print(f"For each prime factor, there are phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {num_choices_div_6} choices.")
    equation_1 = f"{num_choices_div_6}^{len(primes)}"
    print(f"Total combinations where orders divide 6: {equation_1} = {total_div_6}.")

    # 2. Count combinations where the lcm of orders divides 3 (i.e., all orders are 3)
    # For each prime, possible primitive orders are {3}.
    # Number of choices = phi(3) = 2
    num_choices_div_3 = phi(3)
    total_div_3 = num_choices_div_3 ** len(primes)
    print(f"\nStep 2: Subtract combinations where the lcm of orders divides 3 (all orders must be 3).")
    print(f"For each prime factor, there are phi(3) = {num_choices_div_3} choices.")
    equation_2 = f"{num_choices_div_3}^{len(primes)}"
    print(f"Total combinations where orders divide 3: {equation_2} = {total_div_3}.")

    # 3. Count combinations where the lcm of orders divides 2 (i.e., all orders are 2)
    # For each prime, possible primitive orders are {2}.
    # Number of choices = phi(2) = 1
    num_choices_div_2 = phi(2)
    total_div_2 = num_choices_div_2 ** len(primes)
    print(f"\nStep 3: Subtract combinations where the lcm of orders divides 2 (all orders must be 2).")
    print(f"For each prime factor, there are phi(2) = {num_choices_div_2} choices.")
    equation_3 = f"{num_choices_div_2}^{len(primes)}"
    print(f"Total combinations where orders divide 2: {equation_3} = {total_div_2}.")
    
    # Final Calculation
    result = total_div_6 - total_div_3 - total_div_2
    
    print("\nFinal calculation using the Principle of Inclusion-Exclusion:")
    print("The number of characters of order 6 is given by:")
    print(f"(Total combinations with orders dividing 6) - (Total with orders dividing 3) - (Total with orders dividing 2)")
    print(f"Result = {total_div_6} - {total_div_3} - {total_div_2} = {result}")
    
    return result

final_answer = solve()
print(f"\n<<< {final_answer} >>>")