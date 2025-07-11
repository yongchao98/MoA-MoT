import math

def solve():
    """
    Finds the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    order = 6

    # Helper function for prime factorization
    def get_prime_factorization(n):
        factors = {}
        num = 2
        temp = n
        while num * num <= temp:
            while (temp % num) == 0:
                factors[num] = factors.get(num, 0) + 1
                temp //= num
            num += 1
        if temp > 1:
            factors[temp] = factors.get(temp, 0) + 1
        return factors

    print(f"Task: Find the number of primitive Dirichlet characters of conductor d = {d} and order {order}.")
    print("-" * 60)

    # Step 1: Factorize the conductor d and order
    d_factors = get_prime_factorization(d)
    d_primes = sorted(d_factors.keys())
    order_prime_factors = sorted(get_prime_factorization(order).keys())

    print(f"Step 1: Prime factorization of conductor d = {d}.")
    d_factors_str = " * ".join([f"{p}" for p in d_primes])
    print(f"d = {d_factors_str}")
    print(f"This is a square-free conductor.\n")
    print(f"Step 2: Use the Principle of Inclusion-Exclusion for order = {order}.")
    print(f"The distinct prime factors of the order {order} are {order_prime_factors}.")
    print("Let N(k) be the number of primitive characters with conductor d whose order divides k.")
    print("The number of characters of order exactly 6 is N(6) - N(3) - N(2) + N(1).\n")
    
    # Function to calculate N(k)
    def calculate_N_k(k, d_prime_list):
        if k == 0: return 0
        total_prod = 1
        calc_strs = []
        for p in d_prime_list:
            # The number of non-principal characters mod p with order dividing k is gcd(k, p-1) - 1.
            count_p = math.gcd(k, p - 1) - 1
            calc_strs.append(str(count_p))
            total_prod *= count_p
        
        print(f"Calculating N({k}):")
        print(f"  N({k}) = Product over p|d of (gcd({k}, p-1) - 1)")
        print(f"       = {' * '.join(calc_strs)} = {total_prod}")
        return total_prod

    # Step 3: Apply inclusion-exclusion
    total_count = 0
    equation_parts = []
    
    num_order_factors = len(order_prime_factors)
    
    # Iterate through all subsets of the order's prime factors
    for i in range(1 << num_order_factors):
        subset_prod = 1
        num_items_in_subset = 0
        for j in range(num_order_factors):
            if (i >> j) & 1:
                subset_prod *= order_prime_factors[j]
                num_items_in_subset += 1
        
        sign = 1 if num_items_in_subset % 2 == 0 else -1
        k = order // subset_prod
        
        term_val = calculate_N_k(k, d_primes)
        
        total_count += sign * term_val
        
        if sign == 1:
            equation_parts.append(f"+ {term_val}")
        else:
            equation_parts.append(f"- {term_val}")

    # Step 4: Final Result
    print("\n" + "-" * 60)
    print("Final Calculation:")
    final_equation = " ".join(equation_parts).lstrip('+ ')
    print(f"Number of characters = {final_equation}")
    print(f"                     = {total_count}")
    return total_count

# Run the solver to get the final answer.
final_answer = solve()
<<<608>>>