import math

def display_pendulum_period_derivation():
    """
    This function prints the step-by-step derivation and the final answer
    for the period of the pendulum-like system.
    """
    
    # Header
    print("This script presents the derivation of the period of motion for the described system.")
    print("The derivation uses Lagrangian mechanics for small oscillations.")
    
    # Step 1: Equations of Motion
    print("\nStep 1: Equation of Motion")
    print("Using Lagrangian mechanics and applying the small angle approximation, we derive the following")
    print("equation of motion for the angle θ:")
    print("\n  (101/10) * M * R² * d²θ/dt² + 4 * M * g * R * θ = 0\n")
    print("where M is the mass of a disk, R is its radius, and g is the gravitational acceleration.")

    # Step 2: Standard SHM Form
    print("\nStep 2: Re-arrange to Standard Form (θ'' + ω²θ = 0)")
    print("To find the angular frequency (ω), we rearrange the equation:")
    print("\n  d²θ/dt² + ( (4 * g * R) / (101/10 * R²) ) * θ = 0")
    print("  d²θ/dt² + ( (40 * g) / (101 * R) ) * θ = 0\n")

    # Step 3: Identify Angular Frequency and Period
    print("\nStep 3: Identify Angular Frequency (ω) and Period (T)")
    print("From this standard form, we can identify the square of the angular frequency, ω²:")
    print("\n  ω² = (40 * g) / (101 * R)\n")
    print("The period T is defined as 2π / ω. Therefore:")
    print("\n  T = 2π / sqrt( (40 * g) / (101 * R) )")
    
    # Step 4: Final Expression
    print("\n----------------------------------------------------------")
    print("Final Answer: The period of motion T is:")
    print("\n  T = 2 * π * sqrt( (101 * R) / (40 * g) )\n")

    print("Breaking down the final expression as requested:")
    print("Constant Factor = 2")
    print("Constant Factor = π (pi, approx. 3.14159)")
    print("Inside the Square Root:")
    print("  Numerator components: 101 * R")
    print("  Denominator components: 40 * g")
    print("----------------------------------------------------------")

# Execute the function to display the results
display_pendulum_period_derivation()

# The final mathematical expression is T = 2*pi*sqrt(101*R / (40*g))
# Let's extract the core answer string.
final_answer_expression = "2*pi*sqrt(101*R/(40*g))"
# The problem asks for the content of the answer, so the formula is the answer.
print(f'<<<{final_answer_expression}>>>')