def solve_focal_length_exponent():
    """
    This function determines the exponent 'n' in the relationship f ∝ t^n
    for a liquid-mirror telescope spun up by a constant power source.
    """

    print("We are tasked with finding the exponent 'n' in the relationship f ∝ tⁿ.")
    print("Here are the steps to solve the problem:\n")

    # Step 1: Relate focal length (f) to angular velocity (ω)
    print("Step 1: Relate focal length (f) to angular velocity (ω)")
    print("The surface of a rotating fluid forms a paraboloid. Its height 'z' at a radius 'r' is given by:")
    print("z = ω² * r² / (2 * g), where 'g' is the acceleration due to gravity.")
    print("The standard equation for a parabolic mirror with focal length 'f' is:")
    print("z = r² / (4 * f)")
    print("By comparing these two equations, we find the focal length:")
    print("1 / (4 * f) = ω² / (2 * g)  =>  f = g / (2 * ω²)")
    print("This shows that the focal length 'f' is inversely proportional to the square of the angular velocity 'ω'.")
    print("f ∝ ω⁻²\n")

    # Step 2: Relate angular velocity (ω) to time (t) for a constant power source
    print("Step 2: Relate angular velocity (ω) to time (t)")
    print("The system is driven by a constant power source 'P'. In a rotating system, Power P = Torque (τ) * Angular Velocity (ω).")
    print("Torque is also defined as τ = I * α, where 'I' is the moment of inertia and 'α' is the angular acceleration (dω/dt).")
    print("Substituting for τ, we get: P = (I * dω/dt) * ω.")
    print("Since 'P' and 'I' are constant, we can solve this differential equation by separating variables:")
    print("P * dt = I * ω * dω")
    print("We integrate from the initial state (t=0, ω=0) to a general state (t, ω):")
    print("∫[0,t] P dt' = ∫[0,ω] I ω' dω'")
    print("This yields: P * t = (1/2) * I * ω²")
    print("Solving for ω², we find: ω² = (2 * P / I) * t")
    print("This means ω² is directly proportional to time 't'.")
    print("ω² ∝ t\n")

    # Step 3: Combine the relationships to find how 'f' depends on 't'
    print("Step 3: Combine the relationships")
    print("From Step 1, we have the relation: f ∝ ω⁻²")
    print("From Step 2, we have the relation: ω² ∝ t")
    print("By substituting the expression for ω² from Step 2 into the equation from Step 1, we get:")
    print("f ∝ (t)⁻¹\n")

    # Step 4: Identify the exponent 'n' and output the final equation
    print("Step 4: Identify the exponent 'n'")
    print("The problem asks for the value of 'n' in the expression f ∝ tⁿ.")
    n = -1
    print(f"Our derived relationship is f ∝ t^({n}).")
    print(f"Therefore, the value of the exponent n is {n}.")

solve_focal_length_exponent()