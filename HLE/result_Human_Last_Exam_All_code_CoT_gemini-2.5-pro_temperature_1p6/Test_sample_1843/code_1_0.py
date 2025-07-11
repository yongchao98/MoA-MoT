import math

def phi(n):
    """
    Computes Euler's totient function.
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
    return int(result)

def get_prime_factorization_str(n):
    """
    Returns the prime factorization of n as a formatted string.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    
    return " * ".join([f"{p}^{a}" if a > 1 else str(p) for p, a in sorted(factors.items())])

# Problem parameters
N = 36036
order = 6

print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order {order}.")
print("-" * 20)

# Step 1: Prime factorization of N
N_factored_str = get_prime_factorization_str(N)
print(f"Step 1: The prime factorization of N is {N} = {N_factored_str}.")
# N = 2^2 * 3^2 * 7 * 11 * 13
prime_powers = [4, 9, 7, 11, 13]
print(f"The conductor N can be factored into primitive conductors for the prime powers: {prime_powers}.")
print("-" * 20)

print("Step 2: Analyze the number of suitable primitive characters for each prime power conductor.")
print("A character is 'suitable' if it is primitive and its order divides the target order of 6.")
print("-" * 20)

counts = {}

# Conductor 4 (2^2)
pa = 4
print(f"Analyzing conductor p^a = {pa}:")
print("The group of characters mod 4 is cyclic of order phi(4)=2, so it is isomorphic to C_2.")
print("There is 1 primitive character mod 4, and its order is 2.")
print("Since 2 divides 6, this character is suitable.")
num_chars_4 = 1
counts[pa] = num_chars_4
print(f"Number of suitable primitive characters for conductor {pa}: {num_chars_4}.")
print("-" * 10)

# Conductor 9 (3^2)
pa = 9
print(f"Analyzing conductor p^a = {pa}:")
print("The group of characters mod 9 is cyclic of order phi(9)=6, so it is isomorphic to C_6.")
print("Primitive characters mod 9 have orders that do not divide phi(3)=2. So, possible orders are 3 and 6.")
print("Both 3 and 6 divide the target order 6.")
num_order_3 = phi(3)
num_order_6 = phi(6)
num_chars_9 = num_order_3 + num_order_6
counts[pa] = num_chars_9
print(f"Number of primitive characters of order 3: phi(3) = {num_order_3}")
print(f"Number of primitive characters of order 6: phi(6) = {num_order_6}")
print(f"Total number of suitable primitive characters for conductor {pa}: {num_order_3} + {num_order_6} = {num_chars_9}.")
print("-" * 10)

# Conductor 7
pa = 7
print(f"Analyzing conductor p^a = {pa}:")
print("The group of characters mod 7 is cyclic of order phi(7)=6, so it is isomorphic to C_6.")
print("All non-principal characters are primitive. Their orders are 2, 3, and 6.")
print("All these orders (2, 3, 6) divide the target order 6.")
num_order_2 = phi(2)
num_order_3 = phi(3)
num_order_6 = phi(6)
num_chars_7 = num_order_2 + num_order_3 + num_order_6
counts[pa] = num_chars_7
print(f"Number of primitive characters of order 2: phi(2) = {num_order_2}")
print(f"Number of primitive characters of order 3: phi(3) = {num_order_3}")
print(f"Number of primitive characters of order 6: phi(6) = {num_order_6}")
print(f"Total suitable characters for conductor {pa}: {num_order_2} + {num_order_3} + {num_order_6} = {num_chars_7}.")
print("-" * 10)

# Conductor 11
pa = 11
print(f"Analyzing conductor p^a = {pa}:")
print("The group of characters mod 11 is cyclic of order phi(11)=10, so it is isomorphic to C_10.")
print("Primitive characters have orders 2, 5, and 10.")
print("Of these, only order 2 divides the target order 6.")
num_order_2 = phi(2)
num_chars_11 = num_order_2
counts[pa] = num_chars_11
print(f"Number of suitable primitive characters (order 2): phi(2) = {num_order_2}.")
print("-" * 10)

# Conductor 13
pa = 13
print(f"Analyzing conductor p^a = {pa}:")
print("The group of characters mod 13 is cyclic of order phi(12)=12, so it is isomorphic to C_12.")
print("Primitive characters have orders 2, 3, 4, 6, and 12.")
print("Of these, orders 2, 3, and 6 divide the target order 6.")
num_order_2 = phi(2)
num_order_3 = phi(3)
num_order_6 = phi(6)
num_chars_13 = num_order_2 + num_order_3 + num_order_6
counts[pa] = num_chars_13
print(f"Number of primitive characters of order 2: phi(2) = {num_order_2}")
print(f"Number of primitive characters of order 3: phi(3) = {num_order_3}")
print(f"Number of primitive characters of order 6: phi(6) = {num_order_6}")
print(f"Total suitable characters for conductor {pa}: {num_order_2} + {num_order_3} + {num_order_6} = {num_chars_13}.")
print("-" * 20)

print("Step 3: Combine the results.")
print("We form a character by choosing one suitable primitive character for each prime power conductor.")
total = 1
for p_power in sorted(counts.keys()):
    total *= counts[p_power]

print("\nStep 4: Verify the final order.")
print("The order of the combined character is the lcm of the orders of its components.")
print("1. Divisibility by 2: The character for conductor 4 always has order 2, so the lcm is always divisible by 2.")
print("2. Divisibility by 3: The character for conductor 9 always has order 3 or 6, so the lcm is always divisible by 3.")
print("3. Bounded by 6: All component orders are divisors of 6, so their lcm must also divide 6.")
print("Because the lcm is divisible by 2 and 3, it must be divisible by 6. Since it also divides 6, it must be exactly 6.")
print("Therefore, the number of characters of order 6 is the product of the counts of suitable characters.")
print("-" * 20)

print("The final calculation is the product of the number of suitable characters for each factor.")
final_equation = f"Number of characters = {counts[4]} (for 4) * {counts[9]} (for 9) * {counts[7]} (for 7) * {counts[11]} (for 11) * {counts[13]} (for 13)"
print(final_equation)
print(f"                         = {counts[4]} * {counts[9]} * {counts[7]} * {counts[11]} * {counts[13]} = {total}")
