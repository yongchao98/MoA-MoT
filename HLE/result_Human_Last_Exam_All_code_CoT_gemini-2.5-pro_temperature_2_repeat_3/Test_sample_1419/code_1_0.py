import math

def solve_fixed_point_coupling():
    """
    This script derives and prints the leading order expression for the 
    Wilson-Fisher fixed point coupling u* in ϕ⁴ theory near four dimensions.
    """
    
    # Step 1: Define the beta function for ϕ⁴ theory
    print("Step 1: The Beta Function")
    print("In ϕ⁴ theory near d=4 dimensions, the renormalization group flow of the coupling 'u' is governed by the beta function.")
    print("To one-loop order, the beta function is:")
    print("β(u) = -ϵu + B * u²")
    print("Here, ϵ = 4 - d is a small parameter, and B is a constant determined by the theory.\n")

    # Step 2: State the value of the constant B
    print("Step 2: The Coefficient B")
    B_numerator = 3
    B_denominator_str = "16π²"
    print(f"For the single-component ϕ⁴ theory, the one-loop coefficient B = {B_numerator}/{B_denominator_str}.\n")

    # Step 3: Set up the fixed point condition
    print("Step 3: The Fixed Point Condition")
    print("A fixed point, u*, is a value of the coupling where the theory is scale-invariant. This occurs when the beta function is zero:")
    print("β(u*) = 0")
    print(f"So, we must solve the equation: -ϵu* + ({B_numerator}/{B_denominator_str}) * (u*)² = 0\n")

    # Step 4: Solve the equation for the non-trivial fixed point
    print("Step 4: Solving for the Wilson-Fisher Fixed Point")
    print("This equation can be factored as: u* * (-ϵ + ({B_numerator}/{B_denominator_str}) * u*) = 0")
    print("This gives two solutions:")
    print("  - The trivial 'Gaussian' fixed point: u* = 0")
    print("  - The non-trivial 'Wilson-Fisher' fixed point, found by solving:")
    print(f"    -ϵ + ({B_numerator}/{B_denominator_str}) * u* = 0")
    print(f"    ({B_numerator}/{B_denominator_str}) * u* = ϵ")
    print(f"    u* = ϵ / ({B_numerator}/{B_denominator_str})\n")

    # Step 5: Final expression for u*
    print("Step 5: Final Expression")
    print("Rearranging the terms gives the final leading order expression for the fixed point coupling:")
    final_coeff_numerator = 16
    final_coeff_denominator = 3
    
    print(f"u* = ({final_coeff_numerator}π² / {final_coeff_denominator}) * ϵ")

# Execute the function to print the derivation
solve_fixed_point_coupling()