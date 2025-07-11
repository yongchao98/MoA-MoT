import math

def display_fixed_point_expression():
    """
    This function explains and prints the leading order expression for the
    Wilson-Fisher fixed point coupling (u*) in phi^4 theory near d=4 dimensions.
    """
    
    print("Derivation of the Wilson-Fisher Fixed Point Coupling (u*) in phi^4 Theory:")
    print("------------------------------------------------------------------------")
    print("\n1. The theory is analyzed in d = 4 - ε dimensions, where ε is a small parameter.")
    
    print("\n2. The one-loop Renormalization Group (RG) beta function for the dimensionless coupling 'u' is:")
    print("   β(u) = -ε*u + (3 * u²)/(16 * π²)")
    
    print("\n3. A fixed point (u*) is a point where the coupling is scale-invariant, meaning β(u*) = 0:")
    print("   -ε*u* + (3 * (u*)²)/(16 * π²) = 0")
    
    print("\n4. Solving for the non-trivial (Wilson-Fisher) fixed point, we rearrange the equation:")
    print("   (3 * (u*)²)/(16 * π²) = ε*u*")
    print("   3*u* / (16 * π²) = ε")
    
    # Define the numerical coefficients
    numerator = 16
    denominator = 3
    
    print("\n5. The final leading order expression for the fixed point coupling u* is:")
    
    # Final output as requested, printing each number in the equation.
    print(f"\nu* = ({numerator} * π² / {denominator}) * ε")

# Execute the function to display the result
display_fixed_point_expression()
