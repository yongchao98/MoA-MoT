def solve_rod_sliding_problem():
    """
    This function prints the step-by-step derivation for the angle
    at which a tilting rod begins to slide off a table corner.
    """
    print("### Step-by-Step Derivation ###")
    
    print("\n1. Decomposing the Forces")
    print("The problem asks for the angle θ at which the rod begins to slide. This is a static equilibrium problem where we analyze the forces at the moment sliding is about to occur.")
    print("The primary forces acting on the rod at the corner of the table are:")
    print(" - Gravity (Mg) acting vertically downwards at the rod's center of mass.")
    print(" - The normal force (N) exerted by the corner, acting perpendicular to the rod.")
    print(" - The static friction force (f) exerted by the corner, acting parallel to the rod, preventing it from sliding down.")
    
    print("\nWe resolve the force of gravity into components parallel and perpendicular to the rod:")
    print(" - Component perpendicular to the rod: Mg * cos(θ)")
    print(" - Component parallel to the rod: Mg * sin(θ)")

    print("\n2. Applying Equilibrium Conditions")
    print("For the rod to not accelerate through the corner, the normal force must balance the perpendicular component of gravity. So, our first equation is:")
    print("N = Mg * cos(θ)")

    print("\n3. Condition for Sliding")
    print("The rod begins to slide when the component of gravity pulling it down along its length is equal to the maximum possible static friction force.")
    print("The maximum static friction is given by the formula f_max = μ * N, where μ is the coefficient of static friction.")
    print("Therefore, at the point of sliding, the parallel component of gravity equals the maximum friction:")
    print("Mg * sin(θ) = f_max")
    print("Mg * sin(θ) = μ * N")

    print("\n4. Solving for θ")
    print("Now, we substitute our expression for N from step 2 into the equation from step 3:")
    print("Mg * sin(θ) = μ * (Mg * cos(θ))")

    print("\nWe can cancel the term 'Mg' (Mass * gravity) from both sides, as it's non-zero:")
    print("sin(θ) = μ * cos(θ)")

    print("\nTo solve for θ, we divide both sides by cos(θ):")
    print("sin(θ) / cos(θ) = μ")

    print("\nSince sin(θ) / cos(θ) is tan(θ), we get the final equation relating the angle and the coefficient of friction.")
    part1 = "tan(θ)"
    part2 = "μ"
    print(f"The final equation is: {part1} = {part2}")
    
    print("\nThis means the angle θ at which the rod begins to slide is the angle whose tangent is the coefficient of friction μ.")
    
solve_rod_sliding_problem()
