import math

def solve_expected_value():
    """
    Calculates and explains the expected value of T.
    
    The expected value E[T] is given by the series:
    E[T] = sum_{i=1 to inf} i / (3^i * (i+1)!)
    
    This can be rewritten and solved using the Taylor series for e^x.
    """
    
    # The derivation leads to the expression 3 - 2 * e^(1/3)
    e_to_one_third = math.exp(1/3)
    
    # The first sum component: e^(1/3) - 1
    sum1 = e_to_one_third - 1
    
    # The second sum component: 3 * e^(1/3) - 4
    sum2 = 3 * e_to_one_third - 4
    
    # The final expected value E[T] = sum1 - sum2
    expected_value = sum1 - sum2
    
    print("The problem is to find the expected value of T, where T = sum_{i=1 to m} Y_i / 3^i.")
    print("Following the derivation, the expected value E[T] can be expressed as an infinite series:")
    print("E[T] = sum_{i=1 to inf} i / (3^i * (i+1)!)")
    print("\nThis series is evaluated by splitting the term i/((i+1)!) into (1/i!) - (1/(i+1)!).")
    print("This leads to two separate sums:")
    print("Sum 1 = sum_{i=1 to inf} 1 / (3^i * i!)")
    print("Sum 2 = sum_{i=1 to inf} 1 / (3^i * (i+1)!)")
    print(f"\nUsing the Taylor series for e^x, Sum 1 evaluates to e^(1/3) - 1 = {sum1}")
    print(f"Sum 2 evaluates to 3*e^(1/3) - 4 = {sum2}")
    
    print("\nThe final expected value is the difference between these two sums.")
    print(f"E[T] = (e^(1/3) - 1) - (3*e^(1/3) - 4)")
    print(f"E[T] = 3 - 2*e^(1/3)")
    print(f"The numbers in the final equation are 3, 2, and 1/3.")
    
    print(f"\nThe numerical value is approximately: {expected_value}")

solve_expected_value()