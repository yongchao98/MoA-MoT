import math

def find_blowup_condition(x0):
    """
    For the given system of differential equations, this function calculates the
    upper bound for the initial condition y(0) that leads to a blow-up solution,
    given that x(0) = x0 > 1.

    The solution blows up if y(0) is less than a critical value that depends on x(0).
    This critical value is derived from the separatrix of the system's phase space.
    The equation for the separatrix leads to the condition:
    y(0) < sqrt(1 + 2*x(0) - 3*x(0)^(2/3))
    
    Args:
        x0 (float): The initial condition x(0), which must be greater than 1.
    """
    if x0 <= 1:
        print("Error: The initial condition x(0) must be greater than 1.")
        return

    # Constants in the final equation y(0)^2 < c1 + c2*x0 + c3*x0^p
    c1 = 1
    c2 = 2
    c3 = -3
    p = 2.0 / 3.0

    # Calculate the term inside the square root
    y_bound_squared = c1 + c2 * x0 + c3 * (x0**p)

    if y_bound_squared < 0:
        # This case should not happen for x0 > 1 based on the analysis
        print("No real value for the blow-up boundary exists for this x(0).")
    else:
        y_bound = math.sqrt(y_bound_squared)
        
        print("Given the initial condition x(0) = x0 > 1, the solution of the system blows up if y(0) is less than a critical value.")
        print(f"For x(0) = {x0}, this critical value is determined by the equation:")
        print(f"y_crit = sqrt({c1} + {c2}*x0 - {abs(c3)}*x0^({p:.2f}))")
        print(f"y_crit = {y_bound}")
        print("\nTherefore, the solution blows up for y(0) in the interval:")
        print(f"(-infinity, {y_bound})")

if __name__ == '__main__':
    try:
        x0_input = float(input("Enter the initial value for x(0) (must be > 1): "))
        find_blowup_condition(x0_input)
    except ValueError:
        print("Invalid input. Please enter a numerical value.")
