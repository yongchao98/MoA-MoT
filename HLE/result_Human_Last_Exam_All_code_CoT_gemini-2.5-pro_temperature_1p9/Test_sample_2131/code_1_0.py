import math

def calculate_deflection():
    """
    Calculates the membrane's deflection at x=0, which is given by the formula:
    y(0) = (3/4) * (9/4)^(3/5)
    """
    
    # The numbers in the final derived equation for y(0)
    a = 3
    b = 4
    c = 9
    d = 4
    e = 3
    f = 5
    
    # Calculate the value of the expression
    result = (a / b) * math.pow((c / d), (e / f))
    
    # Print the equation with its components and the final result
    print("The final equation for the deflection y(0) is derived as:")
    print(f"y(0) = ({a}/{b}) * ({c}/{d})^({e}/{f})")
    print("\nCalculating this value:")
    print(f"y(0) = {result}")

# Run the calculation
calculate_deflection()
