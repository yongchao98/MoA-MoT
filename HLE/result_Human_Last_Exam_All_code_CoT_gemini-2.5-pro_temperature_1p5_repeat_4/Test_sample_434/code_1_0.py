import math

def solve_blowup_condition():
    """
    This script calculates the critical value for y(0) that determines
    if the solution to the given system of ODEs will blow up,
    assuming x(0) > 1.
    """
    # Set the initial condition for x(t).
    # You can change this value to any number greater than 1.
    x0 = 2.0

    print(f"Analyzing the system of differential equations:")
    print(f"  x'(t) = -3*x(t)*y(t)")
    print(f"  y'(t) = -y(t)^2 - x(t) + 1")
    print(f"for an initial condition x(0) = {x0}, which is greater than 1.")
    print("-" * 50)

    # The solution blows up if y(0) is less than a critical value, yc.
    # This critical value is derived from the stable manifold of the saddle
    # point at (1,0) and is given by the formula:
    # yc = sqrt(2*x(0) + 1 - 3*x(0)^(2/3))

    # Calculate the terms of the equation for clarity
    term1 = 2 * x0
    term2 = 1.0
    term3 = 3 * (x0**(2/3))
    
    # The argument of the square root
    yc_squared = term1 + term2 - term3
    
    if yc_squared < 0:
        # This should not happen for x0 > 1 as proven in the analysis
        print("Error: The argument of the square root is negative.")
        return

    # Calculate the critical value
    yc = math.sqrt(yc_squared)

    print("The solution blows up if y(0) is less than a critical value yc.")
    print(f"The formula for yc is: yc = sqrt(2*x(0) + 1 - 3*x(0)^(2/3))")
    print(f"\nFor x(0) = {x0}, the terms in the equation are:")
    # This printout aims to fulfill the "output each number" requirement
    print(f"  2*x(0)       = {term1:.4f}")
    print(f"  1            = {term2:.4f}")
    print(f"  3*x(0)^(2/3) = {term3:.4f}")
    
    print(f"\nPlugging these into the formula:")
    print(f"yc = sqrt({term1:.4f} + {term2:.4f} - {term3:.4f})")
    print(f"yc = sqrt({yc_squared:.4f})")
    print(f"yc = {yc:.4f}")

    print("-" * 50)
    print(f"CONCLUSION: For x(0) = {x0}, the solution blows up if y(0) < {yc:.4f}")

if __name__ == '__main__':
    solve_blowup_condition()
