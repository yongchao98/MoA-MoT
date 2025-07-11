import math

def solve():
    """
    This function calculates the final value of the integral based on the analytical solution derived above.
    The integral was decomposed into four parts, and each part was evaluated using standard formulas.
    The final answer is of the form A*pi^8 + B*pi^2 + C*pi + D.
    This script prints the coefficients and the final numerical result.
    """

    # Coefficients from the analytical solution
    A = 8/15
    B = 1/3
    C = -1/2
    D = 1

    pi = math.pi
    
    # Calculate the final numerical value
    final_value = A * pi**8 + B * pi**2 + C * pi + D
    
    print("The final result of the integral is given by the equation:")
    print(f"{A} * pi^8 + {B} * pi^2 + {C} * pi + {D}")
    print("\nWhere the numbers in the final equation are:")
    print(f"Coefficient of pi^8: 8/15")
    print(f"Coefficient of pi^2: 1/3")
    print(f"Coefficient of pi: -1/2")
    print(f"Constant term: 1")
    
    print(f"\nThe numerical value is approximately: {final_value}")

solve()