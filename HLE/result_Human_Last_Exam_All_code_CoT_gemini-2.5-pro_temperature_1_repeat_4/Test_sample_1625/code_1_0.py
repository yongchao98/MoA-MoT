import numpy as np

def demonstrate_fourier_series_for_poloidal_angle():
    """
    Demonstrates a Fourier series representation for a function
    dependent on a periodic angle, like the poloidal angle in a torus.
    """
    # In a Fourier series, a function f(theta) is represented as:
    # f(theta) = a0/2 + sum_{n=1 to inf} [ an*cos(n*theta) + bn*sin(n*theta) ]
    # where theta is the poloidal angle (0 to 2*pi).
    
    # Let's define some non-zero Fourier coefficients for our example.
    # These would typically be calculated from experimental data or simulation results.
    # a_n are coefficients for cosine terms, b_n for sine terms.
    # Let's say a1 = 1.5, b2 = -0.8, and a3 = 0.5. All others are zero.
    
    coefficients = {
        'a1': 1.5,
        'b2': -0.8,
        'a3': 0.5
    }
    
    print("In toroidal systems, the poloidal dependence of a quantity is periodic.")
    print("Therefore, it is naturally represented by a Fourier series.")
    print("\nAn example of a function F(θ) representing poloidal variation is:")
    
    equation_parts = []
    # Note: We are ignoring the a0 term for simplicity in this example.
    
    # Process cosine terms (a_n)
    for n in range(1, 5): # Check for a few modes
        coeff_key = f'a{n}'
        if coeff_key in coefficients:
            val = coefficients[coeff_key]
            # Add sign and value
            sign = "+" if val > 0 else "-"
            term = f" {sign} {abs(val)}*cos({n}*θ)"
            equation_parts.append(term)
            
    # Process sine terms (b_n)
    for n in range(1, 5): # Check for a few modes
        coeff_key = f'b{n}'
        if coeff_key in coefficients:
            val = coefficients[coeff_key]
            # Add sign and value
            sign = "+" if val > 0 else "-"
            term = f" {sign} {abs(val)}*sin({n}*θ)"
            equation_parts.append(term)
    
    # Join the parts into a final equation string
    # and clean up the leading " + " if it exists.
    final_equation = "".join(equation_parts).strip()
    if final_equation.startswith('+'):
        final_equation = final_equation[1:].strip()
        
    print(f"F(θ) = {final_equation}")
    print("\nEach term in the equation is a component of the Fourier series expansion.")

demonstrate_fourier_series_for_poloidal_angle()