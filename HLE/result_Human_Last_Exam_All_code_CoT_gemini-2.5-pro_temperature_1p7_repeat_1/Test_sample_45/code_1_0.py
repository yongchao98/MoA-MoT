import sympy

def solve_liquid_mirror_telescope():
    """
    This script derives the relationship between the focal length (f) of a liquid mirror
    and time (t) when spun up by a constant power source. It then computes the exponent
    n in the relationship f ∝ t^n.
    """
    
    print("### Step-by-step Derivation ###")
    
    # --- Step 1: Relate focal length (f) to angular speed (ω) ---
    print("\n--- Step 1: Find the relationship between focal length (f) and angular speed (ω) ---")
    print("For a fluid rotating at a constant angular speed ω, the equilibrium surface shape is a paraboloid.")
    print("The equation for the surface height z at a radius r is given by:")
    print("z(r) = (ω^2 * r^2) / (2 * g), where g is the acceleration due to gravity.")
    print("\nThe standard equation for a parabolic mirror with its vertex at the origin is:")
    print("z(r) = r^2 / (4 * f), where f is the focal length.")
    print("\nBy comparing these two equations, we can relate f and ω:")
    print("(ω^2 * r^2) / (2 * g) = r^2 / (4 * f)")
    print("Solving for f, we find: f = g / (2 * ω^2)")
    print("This shows that the focal length is inversely proportional to the square of the angular speed: f ∝ ω^(-2).")
    
    # --- Step 2: Relate angular speed (ω) to time (t) ---
    print("\n--- Step 2: Find the relationship between angular speed (ω) and time (t) ---")
    print("The system is accelerated from rest by a constant power source, P.")
    print("In rotational mechanics, power (P) is the product of torque (τ) and angular speed (ω): P = τ * ω.")
    print("Torque is also related to the moment of inertia (I) and angular acceleration (α = dω/dt): τ = I * α.")
    print("\nCombining these, we get: P = (I * dω/dt) * ω.")
    print("Since P and I are constant, we can set up a differential equation to solve for ω(t):")
    print("I * ω * dω = P * dt")
    print("We integrate both sides, noting that the rotation starts from rest (ω=0 at t=0):")
    print("∫(I * ω') dω' from 0 to ω = ∫(P) dt' from 0 to t")
    print("The result of the integration is: (1/2) * I * ω^2 = P * t")
    print("\nSolving for ω^2, we get: ω^2 = (2 * P / I) * t.")
    print("This shows that ω^2 is directly proportional to time t: ω^2 ∝ t.")
    
    # --- Step 3: Combine results and find n ---
    print("\n--- Step 3: Combine the relationships to find f as a function of t ---")
    print("From Step 1, we have: f ∝ ω^(-2)")
    print("From Step 2, we have: ω^2 ∝ t")
    print("\nSubstitute the relation for ω^2 into the relation for f:")
    print("f ∝ (ω^2)^(-1)")
    print("f ∝ t^(-1)")
    
    # --- Final Answer ---
    print("\n### Final Result ###")
    print("The relationship is of the form f ∝ t^n.")
    # The final equation is f ∝ t^-1
    n = -1
    print(f"By comparing our result to this form, we find the final equation: f ∝ t^({n})")
    print(f"Therefore, the value of n is {n}.")

solve_liquid_mirror_telescope()

# The final answer is the value of n
final_answer = -1
print(f"\n<<<{final_answer}>>>")