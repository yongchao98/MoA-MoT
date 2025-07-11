import numpy as np

def solve():
    """
    This function calculates the expression for n=1, which corresponds to the least upper bound.
    The formula derived is f(n) = 2 * (n+2) / (n+1).
    The least upper bound (supremum) of f(n) for n>=1 is at n=1.
    """
    n = 1
    
    # E_H = (n+2)/(n+1)
    e_h_num = n + 2
    e_h_den = n + 1
    
    # S_H = n+1
    s_h = n + 1
    
    # E_P = 2/(n+1)
    e_p_num = 2
    e_p_den = n + 1
    
    # S_P = 1
    s_p = 1
    
    # Product = E_P * E_H * S_P * S_H
    # = (2/(n+1)) * ((n+2)/(n+1)) * 1 * (n+1)
    # = 2 * (n+2) / (n+1)
    
    result_num = e_p_num * e_h_num * s_p * s_h
    result_den = e_p_den * e_h_den
    
    # The term (n+1) from S_H cancels one (n+1) in the denominator.
    final_numerator = 2 * (n + 2)
    final_denominator = n + 1

    final_result = final_numerator / final_denominator
    
    print("The problem is to find the least upper bound of the product E_P * E_H * S_P * S_H over all positive integers n.")
    print("The expression for the product as a function of n is: 2 * (n+2) / (n+1).")
    print("This function is decreasing for n >= 1.")
    print("Therefore, the least upper bound is the value at n=1.")
    print(f"For n = 1, the product is 2 * (1+2) / (1+1) = 2 * 3 / 2 = 3.")
    print(f"The equation is E_P * E_H * S_P * S_H = (2 / (1+1)) * ((1+2) / (1+1)) * 1 * (1+1) = {2/(1+1)} * {(1+2)/(1+1)} * {1} * {1+1} = {final_result}")

solve()
<<<3>>>