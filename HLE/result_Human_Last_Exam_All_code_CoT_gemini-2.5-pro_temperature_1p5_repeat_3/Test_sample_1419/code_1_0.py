import math

def calculate_fixed_point_coupling():
    """
    This function calculates and displays the leading order expression for the 
    Wilson-Fisher fixed point coupling in phi^4 theory near d=4 dimensions.
    """
    
    print("Step 1: The Beta Function")
    print("In phi^4 theory near d=4 dimensions (where d = 4 - ε), the one-loop beta function for the dimensionless coupling 'u' is:")
    print("β(u) = -εu + (3 / (16 * π²)) * u² + O(u³)\n")
    
    print("Step 2: The Fixed Point Condition")
    print("A fixed point u* is a value of the coupling where the theory is scale-invariant, meaning β(u*) = 0.")
    print("Setting the leading terms of the beta function to zero gives:")
    print("-εu* + (3 / (16 * π²)) * (u*)² = 0\n")
    
    print("Step 3: Solving for the Fixed Point u*")
    print("This equation has two solutions. One is the trivial Gaussian fixed point, u* = 0.")
    print("The non-trivial solution, known as the Wilson-Fisher fixed point, is found by factoring out u* (for u* ≠ 0):")
    print("-ε + (3 / (16 * π²)) * u* = 0")
    print("Rearranging the terms to solve for u* gives:")
    print("(3 / (16 * π²)) * u* = ε")
    print("u* = ε * (16 * π² / 3)\n")
    
    print("Step 4: Final Expression and Numerical Value")
    print("The leading order expression for the fixed point coupling is u* = (16π²/3)ε.")
    
    # Define the components of the coefficient
    numerator = 16
    denominator = 3
    pi_value = math.pi
    
    # Calculate the coefficient
    coefficient = (numerator * pi_value**2) / denominator
    
    print("\nTo be explicit, the final equation is:")
    print(f"u* = ({numerator} * π² / {denominator}) * ε")
    print(f"   ≈ ({numerator} * {pi_value**2:.5f} / {denominator}) * ε")
    print(f"   ≈ {coefficient:.5f} * ε")

# Execute the function
calculate_fixed_point_coupling()