import math

def solve_fixed_point_coupling():
    """
    This function explains and derives the leading order expression for the
    fixed point coupling u* in phi^4 theory near four dimensions.
    """
    
    print("Derivation of the Wilson-Fisher Fixed Point Coupling u*")
    print("=" * 55)
    
    # Step 1: State the Beta Function
    print("\nStep 1: The one-loop beta function for the dimensionless coupling u in d = 4 - epsilon dimensions is:")
    print("  beta(u) = -epsilon * u + (3 * u^2) / (16 * pi^2)\n")
    
    # Step 2: Apply the Fixed Point Condition
    print("Step 2: A fixed point u* is a value where the coupling does not change with scale, i.e., beta(u*) = 0.")
    print("  -epsilon * u* + (3 * (u*)^2) / (16 * pi^2) = 0\n")
    
    # Step 3: Solve for the fixed point
    print("Step 3: We factor out u* to find the solutions.")
    print("  u* * [-epsilon + (3 * u*) / (16 * pi^2)] = 0\n")
    print("This gives two possible solutions:")
    print("  1. The Gaussian fixed point: u* = 0 (a non-interacting theory).")
    print("  2. The Wilson-Fisher fixed point: -epsilon + (3 * u*) / (16 * pi^2) = 0.\n")
    
    # Step 4: Isolate u* for the non-trivial solution
    print("Step 4: Solving for the non-trivial Wilson-Fisher fixed point u*:")
    print("  (3 * u*) / (16 * pi^2) = epsilon")
    print("  u* = epsilon * (16 * pi^2) / 3\n")
    
    # Final Result
    numerator = 16
    denominator = 3
    
    print("Final Result:")
    print("-----------")
    print("The leading order expression for the fixed point coupling u* is:")
    # The final equation with each number explicitly shown
    final_equation = f"  u* = ({numerator} * pi^2 / {denominator}) * epsilon"
    print(final_equation)
    
    # Also printing the numerical coefficient
    numerical_coefficient = (numerator * math.pi**2) / denominator
    print(f"\nNumerically, this is approximately:")
    print(f"  u* ≈ {numerical_coefficient:.4f} * epsilon")


if __name__ == "__main__":
    solve_fixed_point_coupling()
