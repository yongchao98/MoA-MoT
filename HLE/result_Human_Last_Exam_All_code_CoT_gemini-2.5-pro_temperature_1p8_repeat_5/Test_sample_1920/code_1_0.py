import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    """
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
    Calculates the number of primitive Dirichlet characters of conductor 53599 and order 6.
    """
    d = 53599
    g = 6
    # Prime factorization of d = 53599 is 7 * 13 * 19 * 31
    primes = [7, 13, 19, 31]
    num_primes = len(primes)
    
    print(f"The conductor is d = {d}, and its prime factorization is {primes}.")
    print(f"We are looking for primitive characters of order g = {g}.\n")

    # The possible orders for component characters are divisors of g, excluding 1.
    divisors_of_g = [k for k in range(2, g + 1) if g % k == 0]
    
    # For each prime p_i, the order of a component character must divide gcd(p_i-1, g).
    # In this case, for all p in {7,13,19,31}, gcd(p-1, 6)=6.
    # So possible orders are divisors of 6 > 1, i.e., {2, 3, 6}.
    possible_orders = divisors_of_g
    print(f"For a character to be primitive of conductor {d}, its component characters must have orders from {possible_orders}.\n")
    
    # Number of characters of each possible order
    phi_values = {k: phi(k) for k in possible_orders}
    print("The number of characters for each possible order k is phi(k):")
    for k, v in phi_values.items():
        print(f"  phi({k}) = {v}")
    print()
    
    # Number of primitive characters for each prime factor whose order divides 6
    num_chars_per_prime = sum(phi_values.values())
    print(f"For each prime factor, the number of primitive characters with an order dividing {g} is:")
    print(f"  phi(2) + phi(3) + phi(6) = {phi(2)} + {phi(3)} + {phi(6)} = {num_chars_per_prime}\n")
    
    # Total combinations of characters whose orders divide 6
    total_combinations = num_chars_per_prime ** num_primes
    print(f"The total number of primitive character combinations where each component's order divides {g} is:")
    print(f"  {num_chars_per_prime}^{num_primes} = {total_combinations}\n")
    
    # We subtract combinations where the lcm order is not 6.
    # Case 1: Order is 2. All components must have order 2.
    num_order_2 = phi_values[2] ** num_primes
    print("Subtracting cases where the final order is not 6:")
    print(f"  Number of characters of order 2 = (phi(2))^{num_primes} = {phi_values[2]}^{num_primes} = {num_order_2}")

    # Case 2: Order is 3. All components must have order 3.
    num_order_3 = phi_values[3] ** num_primes
    print(f"  Number of characters of order 3 = (phi(3))^{num_primes} = {phi_values[3]}^{num_primes} = {num_order_3}\n")

    # Final result
    result = total_combinations - num_order_2 - num_order_3

    print("The number of primitive characters of order 6 is the total combinations minus those of order 2 and 3:")
    print(f"Final calculation: {total_combinations} - {num_order_2} - {num_order_3} = {result}")

solve()

print("\nTo summarize the final calculation in one equation:")
p4 = 4
phi2 = phi(2)
phi3 = phi(3)
phi6 = phi(6)
total_per_prime = phi2 + phi3 + phi6
total_comb = total_per_prime**p4
num_o2 = phi2**p4
num_o3 = phi3**p4
final_result = total_comb - num_o2 - num_o3
print(f"({phi2} + {phi3} + {phi6})^{p4} - {phi2}^{p4} - {phi3}^{p4} = {total_per_prime}^{p4} - {num_o2} - {num_o3} = {total_comb} - {num_o2} - {num_o3} = {final_result}")
print("\n<<<608>>>")
