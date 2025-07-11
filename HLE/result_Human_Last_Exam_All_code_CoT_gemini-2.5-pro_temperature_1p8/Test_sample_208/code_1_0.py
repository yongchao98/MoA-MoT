import numpy as np

def solve_circle_packing_14():
    """
    Calculates the radius of a circle tightly packing 14 unit circles.
    The radius R is the largest real root of the polynomial:
    4R^4 - 20R^3 - 34R^2 + 50R + 69 = 0.
    """
    
    # Define the coefficients of the polynomial equation
    c4 = 4
    c3 = -20
    c2 = -34
    c1 = 50
    c0 = 69
    
    coeffs = [c4, c3, c2, c1, c0]
    
    # Print the equation as requested
    print("The radius R is the largest real root of the equation:")
    print(f"{c4}R⁴ + ({c3})R³ + ({c2})R² + {c1}R + {c0} = 0")
    
    # Calculate the roots of the polynomial
    roots = np.roots(coeffs)
    
    # Filter for the real roots (where the imaginary part is effectively zero)
    real_roots = roots[np.isclose(roots.imag, 0)].real
    
    # The radius R is the largest of these real roots
    radius = np.max(real_roots)
    
    # Format the result to 4 significant digits
    # The '.4g' format specifier is used for significant figures.
    formatted_radius = f"{radius:.4g}"
    
    print(f"\nThe radius of the circle is: {formatted_radius}")

# Execute the function to find the answer
solve_circle_packing_14()