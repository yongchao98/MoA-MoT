import math

def count_elements_where_power_is_identity(n, p):
    """
    Calculates the number of elements g in S_n such that g^p = 1.
    For a prime p, this means the cycle decomposition of g can only have 
    cycles of length 1 or p.
    """
    count = 0
    # Iterate through the possible number of cycles of length p, denoted c_p.
    for c_p in range(n // p + 1):
        # The remaining elements must form cycles of length 1 (fixed points).
        c_1 = n - c_p * p
        
        # Use the formula for the size of a conjugacy class.
        # Number of permutations = n! / (c_1! * 1^c_1 * c_p! * p^c_p)
        try:
            numerator = math.factorial(n)
            denominator = math.factorial(c_1) * math.factorial(c_p) * (p ** c_p)
            count += numerator // denominator
        except ValueError:
            # This can happen if c_1 is negative, which it shouldn't be
            # with the loop range, but as a safeguard.
            pass
    return count

def main():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7
    ORDER_C2_GEN = 2
    ORDER_C5_GEN = 5

    # a[n] will store the total number of homomorphisms from G to S_n.
    # t[n] will store the number of transitive homomorphisms from G to S_n.
    a = [0] * (N + 1)
    t = [0] * (N + 1)

    print("--- Starting Calculation ---")
    
    # We iterate from n=1 to N to solve the recurrence relation.
    for n in range(1, N + 1):
        # Calculate a[n], the total number of homomorphisms to S_n.
        # This is the number of elements x with x^2=1 times the number of elements y with y^5=1.
        num_x = count_elements_where_power_is_identity(n, ORDER_C2_GEN)
        num_y = count_elements_where_power_is_identity(n, ORDER_C5_GEN)
        a[n] = num_x * num_y
        
        print(f"\nFor S_{n}:")
        print(f"  Number of elements x with x^{ORDER_C2_GEN}=1: {num_x}")
        print(f"  Number of elements y with y^{ORDER_C5_GEN}=1: {num_y}")
        print(f"  Total homomorphisms a_{n} = {num_x} * {num_y} = {a[n]}")

        # Calculate the sum for intransitive homomorphisms.
        # sum = sum_{k=1 to n-1} C(n-1, k-1) * t_k * a_{n-k}
        intransitive_sum = 0
        for k in range(1, n):
            term = math.comb(n - 1, k - 1) * t[k] * a[n - k]
            intransitive_sum += term
            
        # Calculate t[n], the number of transitive homomorphisms.
        t[n] = a[n] - intransitive_sum
        
        if n > 1:
            print(f"  Number of intransitive homomorphisms = {intransitive_sum}")
        print(f"  Number of transitive homomorphisms t_{n} = {a[n]} - {intransitive_sum} = {t[n]}")

    # Final result calculation.
    # The number of subgroups of index N is t_N / (N-1)!
    factorial_n_minus_1 = math.factorial(N - 1)
    num_subgroups = t[N] // factorial_n_minus_1
    
    print("\n--- Final Result ---")
    print(f"The number of subgroups of index {N} is calculated as t_{N} / ({N}-1)!.")
    print(f"The number of transitive homomorphisms t_{N} is {t[N]}.")
    print(f"The factorial of ({N}-1) is {factorial_n_minus_1}.")
    print("The final equation is:")
    print(f"{t[N]} / {factorial_n_minus_1} = {num_subgroups}")
    
    print(f"<<<{num_subgroups}>>>")

if __name__ == "__main__":
    main()