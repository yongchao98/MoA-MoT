import math

def solve_phi4_fixed_point():
    """
    Calculates and explains the leading order expression for the Wilson-Fisher
    fixed point in phi^4 theory near d=4 dimensions.
    """
    print("This program finds the leading order expression for the non-trivial (Wilson-Fisher) fixed point coupling u*")
    print("in phi^4 theory in d = 4 - epsilon dimensions.\n")

    print("Step 1: The one-loop beta function")
    print("The Renormalization Group (RG) beta function for the dimensionless coupling u, to one-loop order, is:")
    
    # In standard field theory conventions (e.g., Peskin & Schroeder), the coefficient is 3/(16*pi^2).
    # The '3' comes from the combinatorics of the one-loop Feynman diagram.
    B_numerator = 3
    B_denominator_str = "16 * pi^2"
    print(f"   β(u) = -ε*u + ({B_numerator} / {B_denominator_str}) * u^2\n")

    print("Step 2: Finding the fixed point")
    print("A fixed point u* is a value of the coupling where the beta function is zero (β(u*) = 0).")
    print("At this point, the theory is scale-invariant.")
    print("Setting β(u*) = 0 gives the equation:\n")
    print(f"   -ε*u* + ({B_numerator} / {B_denominator_str}) * (u*)^2 = 0\n")

    print("Step 3: Solving the equation")
    print("We can factor the equation to find the possible solutions:")
    print(f"   u* * [-ε + ({B_numerator} / {B_denominator_str}) * u*] = 0\n")
    print("This equation reveals two fixed points:")
    print("   - The Gaussian fixed point: u* = 0 (which is unstable)")
    print("   - The non-trivial Wilson-Fisher fixed point, found by setting the term in the brackets to zero.\n")

    print("Step 4: Deriving the Wilson-Fisher fixed point")
    print("Solving for the non-trivial fixed point:")
    print(f"   -ε + ({B_numerator} / {B_denominator_str}) * u* = 0")
    print(f"   ({B_numerator} / {B_denominator_str}) * u* = ε")
    print(f"   u* = ({B_denominator_str} / {B_numerator}) * ε\n")

    print("Final Answer: The leading order expression for the Wilson-Fisher fixed point coupling is:")
    
    # Display the final equation with all its numerical components
    final_coeff_num_str = "16 * pi^2"
    final_coeff_den = 3
    print(f"   u* = ( {final_coeff_num_str} / {final_coeff_den} ) * ε")

# Run the function to display the derivation
solve_phi4_fixed_point()