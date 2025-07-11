import math

def find_blowup_condition(x0: float):
    """
    For the given system of ODEs, this function determines and prints the
    range of initial values for y(0) that lead to a blow-up, assuming x(0) > 1.

    Args:
        x0 (float): The initial value x(0), which must be greater than 1.
    """
    if x0 <= 1:
        print("This analysis is valid only for x(0) > 1.")
        return

    # The condition for blow-up is that the initial point (x0, y0) must lie
    # below the separatrix curve passing through the saddle point (1, 0).
    # The equation for this separatrix is y^2 = 2x + 1 - 3x^(2/3).
    # Therefore, we require y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).

    # Coefficients and powers in the separatrix equation
    c1 = 2
    c2 = 1
    c3 = -3
    power_numerator = 2
    power_denominator = 3
    power = power_numerator / power_denominator

    print(f"Given the initial condition x(0) = {x0}, where {x0} > 1.")
    print("The solution to the system of differential equations blows up if y(0) is below a critical value.")
    print("This critical value is determined by a separatrix curve.")
    
    # We output each number in the final equation.
    print("\nThe equation for the separatrix is: y^2 = ({} * x) + {} + ({} * x^({}/{}))".format(c1, c2, c3, power_numerator, power_denominator))
    
    # Calculate the numerical value of the upper bound for y(0)
    try:
        # We check that the term inside the square root is positive, which it is for x0 > 1.
        val_inside_sqrt = c1 * x0 + c2 + c3 * (x0**power)
        upper_bound = math.sqrt(val_inside_sqrt)
        
        print(f"\nFor the solution to blow up, the initial condition y(0) must satisfy:")
        print(f"y(0) < sqrt({c1}*x(0) + {c2} - {abs(c3)}*x(0)^({power_numerator}/{power_denominator}))")
        print(f"\nFor x(0) = {x0}, this evaluates to:")
        print(f"y(0) < {upper_bound:.6f}")

    except ValueError:
        # This case should not be reached for x0 > 1.
        print("\nCould not compute the numerical bound because the value inside the square root is negative.")

# Example: Let's test for an initial condition x(0) = 8
# (chosen because its cube root is an integer).
initial_x = 8.0
find_blowup_condition(initial_x)