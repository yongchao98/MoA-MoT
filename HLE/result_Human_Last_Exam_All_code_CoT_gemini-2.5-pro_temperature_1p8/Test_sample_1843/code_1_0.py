import math

def get_prime_factorization(num):
    """Computes the prime factorization of a number."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n):
    """Computes Euler's totient function."""
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
    return int(result)

def get_divisors(n):
    """Gets all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve():
    """
    Finds the number of primitive Dirichlet characters of conductor N and order 6.
    """
    N = 36036
    order = 6

    factors = get_prime_factorization(N)
    print(f"The prime factorization of N = {N} is {factors}.\n")

    total_count = 1
    all_choices = {}

    print("Counting the number of allowed primitive characters for each prime-power factor:")

    # Case for p=2, a=2 (mod 4)
    # The group of characters mod 4 is C_2. The only primitive character has order 2.
    # Since 2 divides 6, this character is allowed.
    count_4 = 1
    all_choices[4] = count_4
    print(f"For p^a = 4: There is 1 primitive character, which has order 2. Since 2 divides 6, it is a valid choice. Choices = {count_4}.")

    # General case for odd primes
    for p, a in factors.items():
        if p == 2:
            continue

        phi_pa = phi(p**a)
        if a > 1:
            phi_prev = phi(p**(a-1))
        else: # a == 1
            phi_prev = phi(1)

        allowed_orders = []
        for k in get_divisors(phi_pa):
            if k % order != 0 and order % k != 0 and k != 1 and order !=1 : # must be a divisor of order=6
                 if order % k != 0 and k % order !=0:
                    continue
            
            # For p^a, primitive characters' orders k must divide phi(p^a).
            # If a > 1, k must NOT divide phi(p^(a-1)).
            # If a == 1, k must be greater than 1.
            is_primitive_order = False
            if a > 1:
                if k % phi_prev != 0:
                    is_primitive_order = True
            else: # a == 1
                if k > 1:
                    is_primitive_order = True
            
            if is_primitive_order and order % k == 0:
                 allowed_orders.append(k)

        count_pa = sum(phi(k) for k in allowed_orders)
        all_choices[p**a] = count_pa
        total_count *= count_pa
        print(f"For p^a = {p**a}: Allowed primitive orders dividing 6 are {allowed_orders}. The number of choices is "
              f"{' + '.join([f'phi({k})' for k in allowed_orders])} = "
              f"{' + '.join([str(phi(k)) for k in allowed_orders])} = {count_pa}.")

    total_count *= count_4

    print("\nTo find the total number of characters, we multiply the number of choices for each factor.")
    print("This is valid because the least common multiple (lcm) of the orders of the component characters is always 6:")
    print("1. The component for conductor 4 must have order 2, so the lcm is divisible by 2.")
    print("2. The component for conductor 9 must have order 3 or 6, so the lcm is divisible by 3.")
    print("3. All component orders are divisors of 6, so their lcm must also divide 6.")
    print("An integer that divides 6 and is divisible by both 2 and 3 must be 6.")

    calculation_str = " * ".join(map(str, all_choices.values()))
    print(f"\nFinal calculation: {calculation_str} = {total_count}")

solve()
<<<100>>>