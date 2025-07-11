import math

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            result -= result // p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result -= result // n
    return int(result)

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def get_primitive_orders(p, a):
    """Returns a list of orders of primitive characters modulo p^a."""
    m = p**a
    if m == 2: return []  # No primitive characters mod 2
    if m == 4: return [2]  # Primitive character mod 4 has order 2
    
    # For p odd, a=1, primitive orders are d | (p-1) with d > 1
    if a == 1:
        return [d for d in get_divisors(p - 1) if d > 1]
    
    # For p odd, a > 1, primitive orders are d | phi(p^a) but d not | phi(p^(a-1))
    phi_m = phi(m)
    phi_prev = phi(p**(a - 1))
    divs_m = get_divisors(phi_m)
    divs_prev = set(get_divisors(phi_prev))
    return [d for d in divs_m if d not in divs_prev]

def calculate_char_count(factors, k):
    """
    Calculates the number of primitive characters with conductor N
    whose order divides k.
    """
    total_count = 1
    print(f"\nCalculating the number of characters with order dividing {k}:")
    for p, a in factors.items():
        m = p**a
        prim_orders = get_primitive_orders(p, a)
        valid_orders = [d for d in prim_orders if k % d == 0]
        count = sum(phi(d) for d in valid_orders)
        print(f"  For conductor {m}, there are {count} choices.")
        total_count *= count
    print(f"Total number of characters with order dividing {k} is {total_count}.")
    return total_count

def solve():
    """Main function to solve the problem."""
    N = 36036
    target_order = 6
    factors = {2: 2, 3: 2, 7: 1, 11: 1, 13: 1}

    print(f"Let N = {N}. We find the number of primitive Dirichlet characters of conductor N and order {target_order}.")
    print(f"The prime factorization of N is 2^2 * 3^2 * 7 * 11 * 13.")
    print("The number of such characters can be found using the Principle of Inclusion-Exclusion.")
    print("Result = (count of chars with order dividing 6) - (dividing 3) - (dividing 2) + (dividing 1)")

    # Calculate A_k for k = 6, 3, 2, 1
    A6 = calculate_char_count(factors, 6)
    A3 = calculate_char_count(factors, 3)
    A2 = calculate_char_count(factors, 2)
    A1 = 0 # Primitive characters have order > 1
    print("\nNumber of characters with order dividing 1 is 0.")

    final_answer = A6 - A3 - A2 + A1
    
    print("\nFinal Calculation:")
    print(f"The number of primitive characters of order 6 is:")
    print(f"{A6} - {A3} - {A2} + {A1} = {final_answer}")

solve()
<<<100>>>