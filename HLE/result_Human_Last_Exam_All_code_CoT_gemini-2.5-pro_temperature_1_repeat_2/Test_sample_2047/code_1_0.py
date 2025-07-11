import numpy as np
import sys

def solve_l(d, lambd):
    """
    Calculates the value of l(d, lambda) based on the derived analytical formula.
    
    The formula is:
    l(d, lambda) = (1/(2*lambda))*(t2**2 - t1**2) + (d-2)*ln(sinc(t1)/sinc(t2))
    where t1 = arccos(sqrt(3/d)) and t2 = arccos(sqrt(2/d)).
    """
    if d < 4:
        raise ValueError("d must be greater than or equal to 4.")
    if lambd < 1:
        raise ValueError("lambda must be greater than or equal to 1.")

    # Step 1: Calculate dot products and angles t1, t2
    try:
        # cos(t1) = x1^T * p0
        cos_t1 = np.sqrt(3.0 / d)
        # cos(t2) = x2^T * p0
        cos_t2 = np.sqrt(2.0 / d)
        
        # t = arccos(cos(t))
        t1 = np.arccos(cos_t1)
        t2 = np.arccos(cos_t2)
    except ValueError as e:
        print(f"Error in mathematical domain for d={d}: {e}")
        return

    # Step 2: Calculate the two main terms of the equation for l(d, lambda)
    
    # First term, related to lambda
    lambda_term_val = (1 / (2 * lambd)) * (t2**2 - t1**2)
    
    # Second term, related to d
    # sinc(x) = sin(x)/x
    sinc_t1 = np.sin(t1) / t1
    sinc_t2 = np.sin(t2) / t2
    
    d_term_val = (d - 2) * np.log(sinc_t1 / sinc_t2)

    # Step 3: Calculate the final result
    result = lambda_term_val + d_term_val
    
    # Output the numbers in the final equation as requested
    print(f"Calculating l(d={d}, lambda={lambd}):")
    print(f"t1 = arccos(sqrt(3/{d})) = {t1:.6f}")
    print(f"t2 = arccos(sqrt(2/{d})) = {t2:.6f}")
    print(f"Final equation: l({d}, {lambd}) = term_lambda + term_d")
    print(f"term_lambda = (1/(2*{lambd}))*({t2:.6f}^2 - {t1:.6f}^2) = {lambda_term_val:.6f}")
    print(f"term_d = ({d}-2)*ln(sinc({t1:.6f})/sinc({t2:.6f})) = {d_term_val:.6f}")
    print(f"l({d}, {lambd}) = {lambda_term_val:.6f} + {d_term_val:.6f} = {result:.6f}")
    
    return result

if __name__ == '__main__':
    # Set default values for d and lambda
    d_val = 4
    lambda_val = 1.0

    # Override with command-line arguments if provided
    if len(sys.argv) == 3:
        try:
            d_val = int(sys.argv[1])
            lambda_val = float(sys.argv[2])
        except ValueError:
            print("Usage: python your_script_name.py <integer d>=4> <float lambda>=1.0>")
            sys.exit(1)
            
    solve_l(d_val, lambda_val)