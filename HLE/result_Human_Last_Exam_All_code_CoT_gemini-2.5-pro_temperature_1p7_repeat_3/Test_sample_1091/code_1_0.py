import math

def calculate_limit():
    """
    This function calculates the analytical limit of n * P(n) as n goes to infinity.
    The derived formula for the limit is (2 * sqrt(3)) / pi.
    """
    
    # The components of the final equation
    two = 2
    sqrt_three = math.sqrt(3)
    pi = math.pi
    
    # The formula for the limit
    limit_value = (two * sqrt_three) / pi
    
    print("The analytical derivation shows that the limit of n * P(n) as n -> infinity is given by the equation:")
    print("Limit = (2 * sqrt(3)) / pi")
    print("-" * 30)
    print("The values of the components are:")
    print(f"The number 2 is: {two}")
    print(f"The square root of 3 is: {sqrt_three}")
    print(f"The value of pi is: {pi}")
    print("-" * 30)
    print("The final result of the limit is:")
    print(limit_value)

if __name__ == '__main__':
    calculate_limit()