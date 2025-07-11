import math

def solve_cubic_root():
    """
    Finds the real root of the equation x^3 - x^2 - 1 = 0 using the bisection method.
    """
    f = lambda x: x**3 - x**2 - 1
    
    # The root is between 1 and 2, as f(1) = -1 and f(2) = 3.
    low = 1.0
    high = 2.0
    
    # Bisection method for finding the root
    for _ in range(100): # 100 iterations is more than enough for high precision
        mid = (low + high) / 2
        if f(mid) < 0:
            low = mid
        else:
            high = mid
            
    return (low + high) / 2

def calculate_minimal_cost_factor():
    """
    Calculates the minimal number of coins per bit of information needed to sort the array.
    """
    print("The optimal strategy relies on solving the cubic equation x^3 - x^2 - 1 = 0.")
    
    # 1. Find the root of the cubic equation
    x_root = solve_cubic_root()
    print(f"The real root of this equation is approximately x = {x_root:.5f}")
    
    print("\nThe cost of sorting is asymptotically given by C * n * log2(n).")
    print("The constant C represents the minimal cost per bit of information.")
    
    # 2. Calculate the constant C = 1 / log2(x)
    print("The constant C is calculated as 1 / log2(x).")
    
    log2_x = math.log2(x_root)
    print(f"First, we calculate log2(x) = log2({x_root:.5f}) = {log2_x:.5f}")
    
    cost_factor = 1 / log2_x
    print(f"Then, we find C = 1 / {log2_x:.5f} = {cost_factor:.5f}")
    
    final_answer = round(cost_factor, 3)
    print(f"\nThus, the minimal cost per bit is approximately {final_answer}.")
    print("The total cost for large n is this factor times the information needed, which is log2(n!) ~ n*log2(n).")

calculate_minimal_cost_factor()
