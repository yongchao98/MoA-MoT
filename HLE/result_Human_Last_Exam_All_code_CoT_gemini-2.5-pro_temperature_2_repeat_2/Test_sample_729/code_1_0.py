import math

def count_power_subgroups_q128():
    """
    Calculates and explains the number of power subgroups in the
    generalized quaternion group of size 128 (Q_128).
    """
    
    print("The generalized quaternion group of size 128 is Q_128.")
    print("Its power subgroups H = Q_128^k = {g^k | g in Q_128} depend on the integer k.\n")

    # Case 1: k is odd
    # For an odd k, Q_128^k = Q_128. This gives one unique power subgroup.
    num_subgroups_k_odd = 1
    print(f"Case 1: k is an odd integer.")
    print(f"For any odd k, the power subgroup Q_128^k is the entire group Q_128.")
    print(f"This case contributes {num_subgroups_k_odd} power subgroup.\n")

    # Case 2: k is even
    # For an even k, Q_128^k = <x^d> where d = gcd(k, 64).
    # Since k is even, d must be an even divisor of 64.
    # We count the number of even divisors of 64.
    n = 64
    divisors = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divisors.add(i)
            divisors.add(n//i)
    
    even_divisors = sorted([d for d in divisors if d % 2 == 0])
    num_subgroups_k_even = len(even_divisors)

    print(f"Case 2: k is an even integer.")
    print(f"For any even k, the power subgroup is of the form <x^d> where d = gcd(k, 64).")
    print(f"This means d must be an even divisor of 64.")
    print(f"The even divisors of 64 are: {even_divisors}")
    print(f"This case contributes {num_subgroups_k_even} distinct power subgroups.\n")
    
    # Total count
    total_power_subgroups = num_subgroups_k_odd + num_subgroups_k_even
    print("The total number of power subgroups is the sum from both cases.")
    print(f"Total Number = (subgroups from odd k) + (subgroups from even k)")
    print(f"Total Number = {num_subgroups_k_odd} + {num_subgroups_k_even} = {total_power_subgroups}")

# Run the calculation
count_power_subgroups_q128()