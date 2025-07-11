import sympy

def solve_telescope_problem():
    """
    This function derives and computes the exponent n in the relation f ∝ t^n
    for a liquid-mirror telescope spun up by a constant power source.
    """

    # --- Step 1: Relate focal length (f) to angular velocity (omega) ---
    # The equation for the surface of a liquid rotating at angular speed ω is:
    # z = ω^2 * r^2 / (2*g)
    # The equation for a parabolic mirror with focal length f is:
    # z = r^2 / (4*f)
    # By comparing these two equations, we can relate f and ω.
    # 1 / (4*f) = ω^2 / (2*g)
    # Solving for f: f = g / (2 * ω^2)
    # This shows that f is proportional to 1/ω^2 or ω^(-2).
    print("Step 1: The focal length (f) is related to the angular velocity (ω) as:")
    print("f ∝ 1 / ω^2  or  f ∝ ω^(-2)")
    print("-" * 30)

    # --- Step 2: Relate angular velocity (ω) to time (t) ---
    # The telescope is spun up by a constant power source, P.
    # The work done (W) by this source over time t is W = P * t.
    # By the work-energy theorem, this work is converted into rotational kinetic energy (K).
    # K = (1/2) * I * ω^2, where I is the moment of inertia (a constant).
    # Since the system starts from rest, the change in kinetic energy is K.
    # So, P * t = (1/2) * I * ω^2.
    # Since P and I are constants, this means ω^2 is directly proportional to t.
    print("Step 2: The angular velocity squared (ω^2) is related to time (t) as:")
    print("ω^2 ∝ t")
    print("-" * 30)

    # --- Step 3: Combine the relationships to find f(t) ---
    # From Step 1: f ∝ 1 / ω^2
    # From Step 2: ω^2 ∝ t
    # Substituting the second relationship into the first gives:
    # f ∝ 1 / t  or  f ∝ t^(-1)
    print("Step 3: Combining the relationships from Step 1 and 2, we find:")
    print("f ∝ 1/t  or  f ∝ t^(-1)")
    print("-" * 30)
    
    # --- Step 4: Determine the exponent n ---
    # We are looking for the exponent n in the relation f ∝ t^n.
    # From our derivation in Step 3, we can see that the exponent is -1.
    n = -1
    print("Step 4: The problem asks for the exponent n in the relation f ∝ t^n.")
    print(f"Our derivation shows that n = {n}.")
    print("-" * 30)
    
    # --- Final Result ---
    print("Final Answer:")
    print("The final relationship between focal length and time is expressed as:")
    print(f"f ∝ t^({n})")


solve_telescope_problem()