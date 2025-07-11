import math

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    # Maximum index we are interested in
    max_index = 7
    
    # Precomputed factorials
    fact = [math.factorial(i) for i in range(max_index + 1)]

    # N2[k] = number of elements x in S_k such that x^2 = 1
    # N5[k] = number of elements y in S_k such that y^5 = 1
    # We can precompute these combinatorial values.
    # N2[k] is the number of involutions (including id) in S_k.
    # N5[k] is the number of elements of order 1 or 5 in S_k.
    N2 = [0, 1, 2, 4, 10, 26, 76, 232]
    N5 = [0, 1, 1, 1, 1, 25, 145, 505]
    
    # H[k] = |Hom(G, S_k)| = N2[k] * N5[k]
    # H[0] is defined as 1 for the formula.
    H = [1] + [N2[k] * N5[k] for k in range(1, max_index + 1)]

    # u[k] will store the number of subgroups of index k
    u = [0] * (max_index + 1)
    
    print("Calculating the number of subgroups u[n] for n = 1 to 7.")
    print("-" * 30)
    
    # Iteratively calculate u[n] for n = 1 to max_index
    for n in range(1, max_index + 1):
        sum_val = 0
        for k in range(1, n):
            sum_val += u[k] * H[n - k] / fact[n - k]
            
        u[n] = round(H[n] / fact[n - 1] - sum_val)
        print(f"Number of subgroups of index {n}, u[{n}] = {u[n]}")

    print("-" * 30)
    print("Detailed calculation for u[7]:")
    
    # Print the formula for u[7]
    n = 7
    formula_str = f"u[{n}] = H[{n}]/fact[{n-1}] - ("
    terms = []
    for k in range(1, n):
        if u[k] > 0: # Only show non-zero terms
            terms.append(f"u[{k}]*H[{n-k}]/fact[{n-k}]")
    formula_str += " + ".join(terms) + ")"
    print(formula_str)

    # Print the formula with numbers plugged in
    val_str = f"u[7] = {H[n]}/{fact[n-1]} - ("
    terms = []
    
    # For common denominator calculation
    numerator_sum = 0
    
    for k in range(1, n):
        if u[k] > 0:
            # Numerator for the sum term
            num_k = u[k] * H[n-k] * fact[n-1] // fact[n-k]
            numerator_sum += num_k
            
            # Formatted string for printing
            terms.append(f"{u[k]}*{H[n-k]}/{fact[n-k]}")

    val_str += " + ".join(terms) + ")"
    print(val_str)
    
    # Perform calculation with a common denominator to maintain precision
    numerator_Hn = H[n]
    denominator = fact[n-1]
    
    final_numerator = numerator_Hn - numerator_sum
    result = final_numerator // denominator
    
    # Print final calculation steps
    print(f"u[7] = ({H[n]} - (", end="")
    num_terms = []
    for k in range(1,n):
        if u[k]>0:
             num_terms.append(f"{u[k] * H[n - k] * fact[n - 1] // fact[n - k]}")
    print(" + ".join(num_terms) + f")) / {denominator}")
    print(f"u[7] = ({H[n]} - {numerator_sum}) / {denominator}")
    print(f"u[7] = {final_numerator} / {denominator}")
    print(f"The number of subgroups of index 7 is {result}.")
    
solve()
