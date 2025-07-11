import math

def solve_fixed_point_coupling():
    """
    This function derives and prints the leading order expression for the 
    Wilson-Fisher fixed point coupling u* in phi^4 theory near four dimensions.
    """
    
    print("This script derives the leading order expression for the fixed point coupling u* in ϕ⁴ theory.")
    print("The method used is the Renormalization Group in d = 4 - ϵ dimensions.")
    print("-" * 70)

    # Step 1: The beta function
    print("\nStep 1: State the one-loop beta function")
    print("The beta function, β(u), describes the running of the dimensionless coupling 'u' with the energy scale.")
    print("In d = 4 - ϵ dimensions, the one-loop beta function for ϕ⁴ theory is:")
    print("β(u) = -ϵ*u + (3 / (16 * π²)) * u²\n")

    # Step 2: The fixed point condition
    print("Step 2: Find the fixed point u* by setting the beta function to zero")
    print("A fixed point is a scale-invariant point where β(u*) = 0. This gives the equation:")
    print("-ϵ*u* + (3 / (16 * π²)) * (u*)² = 0\n")

    # Step 3: Solve for the non-trivial fixed point
    print("Step 3: Solve the equation for the non-trivial fixed point (the Wilson-Fisher fixed point)")
    print("We can factor out u* to get: u* * [-ϵ + (3 / (16 * π²)) * u*] = 0")
    print("Ignoring the trivial solution u* = 0, we solve the expression in the brackets:")
    print("-ϵ + (3 / (16 * π²)) * u* = 0")
    print("=> (3 / (16 * π²)) * u* = ϵ")
    print("=> u* = ϵ / (3 / (16 * π²))\n")

    # Final Result
    print("Final Result:")
    print("The leading order expression for the fixed point coupling u* is:")
    
    # Define the components of the expression to be printed
    coefficient_numerator = 16
    pi_squared_symbol = "π²"
    coefficient_denominator = 3
    epsilon_symbol = "ϵ"
    
    print(f"u* = ({coefficient_numerator} * {pi_squared_symbol} / {coefficient_denominator}) * {epsilon_symbol}")

if __name__ == "__main__":
    solve_fixed_point_coupling()