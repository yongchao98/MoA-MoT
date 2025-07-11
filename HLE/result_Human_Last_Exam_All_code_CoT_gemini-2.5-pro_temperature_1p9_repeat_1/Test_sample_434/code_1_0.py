import sys
import numpy as np

def calculate_y0_threshold(x0):
    """
    Calculates the critical threshold for y(0) for a given x(0) > 1 based on the separatrix equation.
    The separatrix that dictates the behavior is the stable manifold of the saddle point at (1,0).
    Its equation is y^2 = 2*x + 1 - 3*x^(2/3).
    Solutions blow up for initial conditions (x0, y0) that lie "below" this curve.
    
    Args:
        x0 (float): The initial value for x(t), must be greater than 1.

    Returns:
        A tuple containing the threshold for y0 and the coefficients from the separatrix equation.
    """
    if x0 <= 1:
        raise ValueError("x(0) must be greater than 1 for this analysis.")

    # These are the numbers in the final equation for the separatrix
    coeff_x = 2
    constant = 1
    coeff_x_pow = 3
    power_numerator = 2
    power_denominator = 3
    
    power = power_numerator / power_denominator

    # The equation for the square of the separatrix at x0
    y_sep_sq = coeff_x * x0 + constant - coeff_x_pow * x0**power
    
    if y_sep_sq < 0:
        # This case should not be reached for x0 > 1 since the separatrix function
        # y_sep_sq(x) is 0 at x=1 and its derivative is positive for x>1.
        # This is a safeguard. A negative result means the separatrix isn't real,
        # which shouldn't happen.
        # We can effectively say any y(0) will lead to blow up in this unlikely scenario.
        return float('-inf'), coeff_x, constant, coeff_x_pow, power_numerator, power_denominator


    y_threshold = np.sqrt(y_sep_sq)
    return y_threshold, coeff_x, constant, coeff_x_pow, power_numerator, power_denominator

def main():
    """
    Main function to get user input for x(0) and print the condition for blow-up.
    """
    print("This script determines the condition on y(0) for a solution to blow up, given x(0) > 1.")
    
    try:
        # Check for command-line arguments, otherwise use an example.
        if len(sys.argv) > 1:
            x0 = float(sys.argv[1])
        else:
            x0 = 8.0 
            print(f"No value for x(0) provided via command line. Using an example: x(0) = {x0}\n")
            
        y_thresh, c_x, c, c_xp, p_n, p_d = calculate_y0_threshold(x0)
        
        print("The system blows up if the initial condition (x(0), y(0)) is below the separatrix.")
        # Printing the numbers in the final equation as requested
        print(f"The separatrix is defined by the equation: y^2 = {c_x}*x + {c} - {c_xp}*x^({p_n}/{p_d})")
        
        print(f"\nFor the given initial condition x(0) = {x0}, the solution blows up if:")
        print(f"y(0) < {y_thresh:.6f}")

    except ValueError as e:
        print(f"Error: {e}")
        print("Please provide a numeric value for x(0) greater than 1.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()