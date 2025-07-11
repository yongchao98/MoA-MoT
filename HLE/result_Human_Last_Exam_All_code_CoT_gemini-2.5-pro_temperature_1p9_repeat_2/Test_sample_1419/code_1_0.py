import math

def solve_fixed_point_coupling():
    """
    Calculates and explains the derivation of the leading order expression for the
    Wilson-Fisher fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """

    print("Derivation of the Wilson-Fisher Fixed Point Coupling u*")
    print("======================================================")

    # Step 1: Explain the beta function
    print("\nStep 1: The Beta Function β(u)")
    print("In the renormalization group analysis of ϕ⁴ theory in d = 4 - ϵ dimensions,")
    print("the beta function describes the running of the dimensionless coupling constant 'u'.")
    print("To one-loop order, the beta function is:")
    print("\n  β(u) = -ϵ*u + (3 / (16 * π²)) * u²\n")
    print("The first term, -ϵ*u, arises from the classical scaling of the coupling.")
    print("The second term, proportional to u², is the one-loop quantum correction.")

    # Step 2: Explain the fixed point condition
    print("\nStep 2: The Fixed Point Condition")
    print("A fixed point u* is a value of the coupling where the theory becomes scale-invariant.")
    print("This occurs when the coupling stops 'running', which means the beta function is zero:")
    print("\n  β(u*) = 0\n")

    # Step 3: Solve the equation
    print("\nStep 3: Solving for u*")
    print("We set the beta function expression to zero to find u*:")
    print("  -ϵ*u* + (3 / (16 * π²)) * (u*)² = 0")
    print("\nFactoring out u*, we get:")
    print("  u* * [ -ϵ + (3 / (16 * π²)) * u* ] = 0")
    print("\nThis equation has two solutions:")
    print("1. u* = 0 (the trivial Gaussian fixed point)")
    print("2. -ϵ + (3 / (16 * π²)) * u* = 0 (the non-trivial Wilson-Fisher fixed point)")
    print("\nSolving the second equation for u*:")
    print("  (3 / (16 * π²)) * u* = ϵ")
    print("\nRearranging the terms gives the final expression for the Wilson-Fisher fixed point.")
    
    # Final Answer Section
    numerator = 16
    denominator = 3
    
    print("\n-----------------------------------------------------------")
    print("The leading order expression for the fixed point coupling is:")
    print(f"\n  u* = ({numerator} * π² / {denominator}) * ϵ\n")
    print("-----------------------------------------------------------")

# Run the function to display the derivation
solve_fixed_point_coupling()