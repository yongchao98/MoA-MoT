import sympy

def calculate_theoretical_stress():
    """
    Calculates and explains the theoretical stress at the tip of a sharp wedge
    in a plate under uniaxial tension using symbolic mathematics.
    """
    # Define symbolic variables for the equation.
    # sigma_y: nominal applied stress on the ends of the plate.
    # d: depth of the notch (half the base of the triangular wedge).
    # rho: radius of curvature at the tip of the wedge.
    sigma_y, d, rho = sympy.symbols('σ_y d ρ', positive=True, real=True)

    # Step 1: Define the equation for the maximum stress at the notch.
    # The equation is σ_max = Kt * σ_y, where Kt is the stress concentration factor.
    # For a notch, Kt is approximated as: 1 + 2 * sqrt(d / ρ).
    
    # We can write the full equation for the maximum stress (σ_max).
    # Each part of the equation is represented by a symbol.
    one = 1
    two = 2
    kt_expression = one + two * sympy.sqrt(d / rho)
    sigma_max_eq = kt_expression * sigma_y

    print("The theoretical stress at the tip of the wedge is calculated based on the theory of stress concentration.")
    print("\n--- Step 1: The Governing Equation ---")
    print("The maximum stress (σ_max) at the tip of a notch is:")
    print("σ_max = Kt * σ_y")
    print("where 'σ_y' is the nominal applied stress and 'Kt' is the stress concentration factor.")
    print("\nFor a notch, the factor 'Kt' is approximately:")
    print(f"Kt = {one} + {two} * sqrt(d / ρ)")
    print("where 'd' is the notch depth and 'ρ' is the tip radius.")
    
    print("\nCombining these gives the final equation for stress at the tip:")
    print(f"σ_max = ({one} + {two} * sqrt(d / ρ)) * σ_y")
    print("-" * 60)

    # Step 2: Apply the condition for a theoretically sharp wedge tip.
    # A perfectly sharp corner implies its radius of curvature, ρ, is zero.
    print("\n--- Step 2: Applying the Geometry Condition ---")
    print("For a theoretically sharp corner like the wedge tip (point A), the radius of curvature 'ρ' approaches 0.")
    print("ρ → 0")
    print("-" * 60)

    # Step 3: Calculate the limit to find the theoretical stress.
    # We use sympy's limit function to evaluate the equation as ρ -> 0.
    theoretical_stress_value = sympy.limit(sigma_max_eq, rho, 0, dir='+')

    print("\n--- Step 3: Calculating the Theoretical Stress ---")
    print("We find the value of σ_max by taking the limit of the equation as ρ → 0:")
    print(f"Theoretical Stress = lim(ρ→0) [({one} + {two} * sqrt(d / ρ)) * σ_y]")
    print("\nAs 'ρ' becomes infinitesimally small, the term 'd/ρ' becomes infinitely large,")
    print("causing the entire expression to approach infinity.")
    print("\nThe final theoretical stress value is:")
    print(f"Theoretical Stress = {theoretical_stress_value}")
    
    print("\n\nNote: In reality, stress cannot be infinite. The material would yield or fracture before reaching this theoretical value. However, the theoretical value based on linear elasticity is indeed infinity.")

calculate_theoretical_stress()