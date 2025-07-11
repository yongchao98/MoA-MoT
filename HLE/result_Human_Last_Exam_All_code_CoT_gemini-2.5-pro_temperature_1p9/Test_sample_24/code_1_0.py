import math

# Memoization dictionaries for our recursive functions
A_memo = {}
B_memo = {}

def A(n):
    """
    Calculates the number of elements of order dividing 2 (involutions) in S_n.
    This follows the recurrence A(n) = A(n-1) + (n-1)*A(n-2).
    """
    if n in A_memo:
        return A_memo[n]
    if n == 0:
        return 1
    if n == 1:
        return 1
    result = A(n-1) + (n-1) * A(n-2)
    A_memo[n] = result
    return result

def B(n):
    """
    Calculates the number of elements of order dividing 5 in S_n.
    This follows the recurrence B(n) = B(n-1) + (n-1)(n-2)(n-3)(n-4)*B(n-5).
    """
    if n in B_memo:
        return B_memo[n]
    if n < 5:
        return 1
    
    prod_term = (n-1)*(n-2)*(n-3)*(n-4)
    result = B(n-1) + prod_term * B(n-5)
    B_memo[n] = result
    return result

def combinations(n, k):
    """Helper function for combinations"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n-k))

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    n_max = 7
    N = [0] * (n_max + 1)
    t = [0] * (n_max + 1)
    u = [0] * (n_max + 1)
    
    # Base case
    N[0] = 1

    print("Step 1: Calculate N_k = |Hom(G, S_k)| for k=1 to 7")
    for k in range(1, n_max + 1):
        N[k] = A(k) * B(k)
        print(f"N_{k} = A({k})*B({k}) = {A(k)} * {B(k)} = {N[k]}")

    print("\nStep 2: Calculate t_k (number of transitive homomorphisms) and u_k (number of subgroups)")
    for n in range(1, n_max + 1):
        sum_term = 0
        for k in range(1, n):
            term = combinations(n-1, k-1) * t[k] * N[n-k]
            sum_term += term
        
        t[n] = N[n] - sum_term
        u[n] = t[n] // math.factorial(n-1)
        print(f"t_{n} = {t[n]}, u_{n} = {t[n]} / {n-1}! = {u[n]}")

    print("\nStep 3: Final calculation for the number of subgroups of index 7")
    sum_terms_str = []
    sum_terms_val = 0
    
    # For a more detailed printout of the t_7 calculation
    for k in range(1, 7):
        comb = combinations(6, k - 1)
        term_val = comb * t[k] * N[7 - k]
        sum_terms_val += term_val
        if t[k] != 0:
            sum_terms_str.append(f"{comb} * {t[k]} * {N[7 - k]}")

    t7_sum_calc = N[7] - sum_terms_val
    
    print("\nThe number of transitive homomorphisms t_7 is calculated using the recurrence:")
    print(f"t_7 = N_7 - ( sum_k=1^6 C(6, k-1) * t_k * N_(7-k) )")
    print(f"t_7 = {N[7]} - ( {' + '.join(sum_terms_str).replace('* 0 *','*...*')} )")

    # Manually compute the sum to show the step
    s1 = combinations(6,0)*t[1]*N[6]
    s2 = combinations(6,1)*t[2]*N[5]
    s3 = combinations(6,2)*t[3]*N[4] # this is 0
    s4 = combinations(6,3)*t[4]*N[3] # this is 0
    s5 = combinations(6,4)*t[5]*N[2]
    s6 = combinations(6,5)*t[6]*N[1]
    total_sum = s1+s2+s3+s4+s5+s6
    print(f"t_7 = {N[7]} - ( {s1} + {s2} + 0 + 0 + {s5} + {s6} )")
    print(f"t_7 = {N[7]} - {total_sum}")
    print(f"t_7 = {t[7]}")

    print("\nThe number of subgroups of index 7 is u_7 = t_7 / (7-1)!")
    fact6 = math.factorial(6)
    result = t[7] // fact6
    print(f"u_7 = {t[7]} / {fact6} = {result}")
    
    return result

# Run the calculation and print the final answer in the desired format
final_answer = solve()
print(f"\nFinal Answer: {final_answer}")
print(f"<<<{final_answer}>>>")
