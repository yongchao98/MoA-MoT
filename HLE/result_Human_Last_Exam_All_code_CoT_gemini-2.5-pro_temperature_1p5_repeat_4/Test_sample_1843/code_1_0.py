import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
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
    # A simple implementation is sufficient for the numbers in this problem
    if n == 2: return 1
    if n == 3: return 2
    if n == 4: return 2
    if n == 5: return 4
    if n == 6: return 2
    if n == 7: return 6
    if n == 9: return 6
    if n == 11: return 10
    if n == 13: return 12
    # Fallback for other numbers
    factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return result

# Main logic
N = 36036
print(f"We want to find the number of primitive Dirichlet characters of conductor N = {N} and order 6.")
print("\nStep 1: Find the prime factorization of N.")
factors = get_prime_factorization(N)
factor_str = " * ".join([f"{p}^{a}" for p, a in factors.items()])
print(f"The prime factorization of {N} is: {factor_str}")
print("A character mod N is primitive if and only if it is a product of primitive characters for each prime power factor.")
print("\nStep 2: For each prime power factor, count the number of primitive characters whose order divides 6.")

# Factor q = 4 = 2^2
q4 = 4
print(f"\nFor q = {q4}:")
print(f"The group of characters mod {q4} has order phi({q4}) = {phi(q4)}.")
print(f"There is 1 primitive character modulo {q4}, and its order is 2.")
print("Since 2 divides 6, this is a valid choice.")
count_4 = 1
print(f"Number of choices for characters mod {q4}: {count_4}")

# Factor q = 9 = 3^2
q9 = 9
print(f"\nFor q = {q9}:")
print(f"The group of characters mod {q9} is cyclic of order phi({q9}) = {phi(q9)}.")
print("The primitive characters mod 9 have orders that divide 6 but not phi(3)=2. These orders are 3 and 6.")
count_9 = phi(3) + phi(6)
print(f"Number of characters of order 3 is phi(3)={phi(3)}. Number of characters of order 6 is phi(6)={phi(6)}.")
print(f"Number of choices for characters mod {q9}: {phi(3)} + {phi(6)} = {count_9}")

# Factor q = 7
q7 = 7
print(f"\nFor q = {q7}:")
print(f"The group of characters mod {q7} is cyclic of order phi({q7}) = {phi(q7)}.")
print("For a prime modulus, all non-principal characters are primitive. Their orders are the divisors of 6 (>1): 2, 3, 6.")
count_7 = phi(2) + phi(3) + phi(6)
print(f"Number of characters of order 2 is phi(2)={phi(2)}. Order 3: phi(3)={phi(3)}. Order 6: phi(6)={phi(6)}.")
print(f"Number of choices for characters mod {q7}: {phi(2)} + {phi(3)} + {phi(6)} = {count_7}")

# Factor q = 11
q11 = 11
print(f"\nFor q = {q11}:")
print(f"The group of characters mod {q11} has order phi({q11}) = {phi(q11)}.")
print("The primitive characters have orders dividing 10 (>1), which are 2, 5, 10. Only order 2 divides 6.")
count_11 = phi(2)
print(f"Number of primitive characters of order 2 is phi(2) = {count_11}.")
print(f"Number of choices for characters mod {q11}: {count_11}")

# Factor q = 13
q13 = 13
print(f"\nFor q = {q13}:")
print(f"The group of characters mod {q13} has order phi({q13}) = {phi(q13)}.")
print("The primitive characters have orders dividing 12 (>1). Those dividing 6 are 2, 3, 6.")
count_13 = phi(2) + phi(3) + phi(6)
print(f"Their numbers are phi(2)={phi(2)}, phi(3)={phi(3)}, and phi(6)={phi(6)}.")
print(f"Number of choices for characters mod {q13}: {phi(2)} + {phi(3)} + {phi(6)} = {count_13}")

print("\nStep 3: Check the lcm condition and calculate the total number.")
print("The order of the final character is the lcm of the orders of the components (o_4, o_9, o_7, o_11, o_13).")
print("1. lcm is a multiple of 2: The character mod 4 must have order 2. So, lcm is divisible by 2.")
print("2. lcm is a multiple of 3: The character mod 9 must have order 3 or 6. So, lcm is divisible by 3.")
print("Since the lcm of the orders is always divisible by 2 and 3, and all component orders divide 6, the lcm must be exactly 6.")
print("Thus, the total number of such characters is the product of the number of choices for each component.")

final_result = count_4 * count_9 * count_7 * count_11 * count_13
print("\nFinal calculation:")
print(f"Total number = (choices for mod 4) * (choices for mod 9) * (choices for mod 7) * (choices for mod 11) * (choices for mod 13)")
print(f"Total number = {count_4} * {count_9} * {count_7} * {count_11} * {count_13} = {final_result}")