import math

def solve():
    """
    Calculates the minimum q for Alice to have a 95% win probability
    and returns floor(100*q).
    """
    # Alice's target winning probability
    x = 0.95

    # The formula to find q is derived as:
    # q = [ (1 - (1 - x)^(1/3)) / x^3 ]^(1/4)
    # The code below calculates this value step-by-step.

    # Calculate the term (1-x)
    val_1_minus_x = 1 - x
    
    # Calculate (1-x)^(1/3)
    val_cbrt = math.pow(val_1_minus_x, 1/3)
    
    # Calculate x^3
    val_x_cubed = math.pow(x, 3)
    
    # Calculate the numerator of the fraction inside the root
    numerator = 1 - val_cbrt
    
    # The denominator is x^3
    denominator = val_x_cubed
    
    # Calculate q^4
    q_pow_4 = numerator / denominator
    
    # Calculate q (which is q_0)
    q0 = math.pow(q_pow_4, 1/4)
    
    # Calculate the final result as per the problem statement
    result = math.floor(100 * q0)
    
    print("The equation for q is: q = [ (1 - (1 - x)^(1/3)) / x^3 ]^(1/4)")
    print(f"1. With x = {x}, the term (1 - x) is: {val_1_minus_x}")
    print(f"2. The cube root of that, (1 - x)^(1/3), is: {val_cbrt}")
    print(f"3. The denominator term, x^3, is: {val_x_cubed}")
    print(f"4. The full fraction, q^4, is ({numerator}) / ({denominator}) = {q_pow_4}")
    print(f"5. The value of q is the 4th root of that: {q0}")
    print(f"6. The value of 100*q is: {100 * q0}")
    print(f"7. The floor of 100*q is: {result}")

solve()