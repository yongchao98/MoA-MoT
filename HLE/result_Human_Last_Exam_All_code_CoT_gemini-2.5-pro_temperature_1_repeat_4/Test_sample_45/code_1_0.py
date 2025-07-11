def solve_focal_length_exponent():
    """
    This script derives the relationship between the focal length (f) of a
    liquid-mirror telescope and time (t) when driven by a constant power source.
    It then computes the exponent 'n' in the relation f is proportional to t^n.
    """

    print("Step 1: Find the focal length (f) in terms of angular velocity (ω).")
    print("The surface of a liquid rotating at angular velocity ω forms a paraboloid with the equation z = ω²r² / (2g).")
    print("The equation for a standard parabolic mirror is r² = 4fz.")
    print("By comparing these two equations, we can find the focal length f:")
    print("4f = 2g / ω²")
    print("f = g / (2ω²)")
    print("-" * 50)

    print("Step 2: Find the relationship between angular velocity (ω) and time (t).")
    print("The system is driven by a constant power source (P). Power is the rate of change of kinetic energy (K): P = dK/dt.")
    print("The rotational kinetic energy is K = (1/2)Iω², where I is the moment of inertia.")
    print("So, P = d/dt[(1/2)Iω²]. Since I is constant, P = (1/2)I * d/dt(ω²).")
    print("We integrate this with respect to time, starting from rest (ω=0 at t=0):")
    print("∫ P dt = ∫ (1/2)I d(ω²)")
    print("Pt = (1/2)Iω²")
    print("Solving for ω², we get: ω² = 2Pt / I")
    print("-" * 50)

    print("Step 3: Combine the equations to find f as a function of t.")
    print("We substitute the expression for ω² from Step 2 into the equation for f from Step 1.")
    print("f = g / (2 * [2Pt / I])")
    print("f = gI / (4Pt)")
    print("f = (gI / 4P) * (1/t)")
    print("f = (gI / 4P) * t⁻¹")
    print("-" * 50)

    print("Step 4: Determine the exponent n.")
    print("We are looking for n where f is proportional to tⁿ.")
    print("From our derivation, f is proportional to t⁻¹.")
    final_exponent = -1
    print(f"The final relationship is f ∝ t^({final_exponent})")
    print(f"Therefore, the value of n is {final_exponent}.")


solve_focal_length_exponent()
<<< -1 >>>