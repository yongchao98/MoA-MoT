def solve_raindrop_problem():
    """
    This script derives the formula for the time it takes a raindrop to fall
    while accumulating mass.
    """

    print("Step 1: Setting up the equation of motion.")
    print("Let 'm' be the mass of the raindrop, 'v' its velocity, and 'g' the acceleration due to gravity.")
    print("The raindrop has density 'rho' (ρ) and falls through micro-droplets of density 'Rho' (Ρ).")
    print("Newton's second law for a variable-mass system where the added mass is initially stationary is F_ext = d(mv)/dt.")
    print("The external forces are gravity (F_g = mg) and buoyancy (F_b = -Ρ * V * g).")
    print("The volume of the drop V = m/ρ. So, F_b = -Ρ * (m/ρ) * g.")
    print("The net external force is F_net = mg - (Ρ/ρ)mg = mg(1 - Ρ/ρ).")
    print("So, the equation of motion is: d(mv)/dt = mg(1 - Ρ/ρ)\n")

    print("Step 2: Modeling the mass accumulation.")
    print("Let 'r' be the radius of the spherical raindrop. Its cross-sectional area is A = πr².")
    print("As it falls at velocity 'v', it sweeps a volume 'Av' per second, accumulating mass at a rate dm/dt.")
    print("dm/dt = Ρ * A * v = Ρ * π * r² * v\n")

    print("Step 3: Relating mass, radius, and distance.")
    print("The mass 'm' is related to the radius 'r' by m = ρ * Volume = ρ * (4/3)πr³.")
    print("Taking the time derivative: dm/dt = ρ * 4πr² * (dr/dt).")
    print("Equating the two expressions for dm/dt: Ρπr²v = 4ρπr²(dr/dt).")
    print("This simplifies to v = (4ρ/Ρ) * (dr/dt).")
    print("Since v = dy/dt (where 'y' is the distance fallen), we can integrate to relate 'y' and 'r'.")
    print("Assuming the raindrop starts from zero size and velocity (r(0)=0, y(0)=0), we get y = (4ρ/Ρ)r.\n")

    print("Step 4: Solving the equation of motion.")
    print("Let's substitute everything into d(mv)/dt = mg(1 - Ρ/ρ).")
    print("The equation can be transformed into one for acceleration 'a' as a function of distance 'y':")
    print("a = g(1 - Ρ/ρ) - (3/y)v²")
    print("We test for a constant acceleration solution, a = C. For constant acceleration, v² = 2ay.")
    print("a = g(1 - Ρ/ρ) - (3/y)(2ay) = g(1 - Ρ/ρ) - 6a")
    print("7a = g(1 - Ρ/ρ)")
    print("a = (g/7) * (1 - Ρ/ρ) = (g/7) * (ρ - Ρ)/ρ")
    print("So, the raindrop falls with a constant acceleration.\n")

    print("Step 5: Calculating the time to fall.")
    print("The time 't' to fall a height 'h' under constant acceleration 'a' from rest is given by h = (1/2)at².")
    print("So, t = sqrt(2h/a).")
    print("Substituting our expression for 'a':")
    print("t = sqrt(2h / [(g/7) * (ρ - Ρ)/ρ])")
    print("t = sqrt((14 * h * ρ) / (g * (ρ - Ρ)))\n")
    
    print("Step 6: Substituting the given height and simplifying.")
    print("The problem states the initial height is h = 350g meters.")
    initial_h_factor = 350
    print(f"h = {initial_h_factor} * g")
    print("t = sqrt((14 * (350 * g) * ρ) / (g * (ρ - Ρ)))")
    print("The 'g' in the numerator and denominator cancels out.")
    intermediate_factor = 14 * initial_h_factor
    print(f"t = sqrt(({intermediate_factor} * ρ) / (ρ - Ρ))")
    final_factor = int(intermediate_factor**0.5)
    print(f"t = {final_factor} * sqrt(ρ / (ρ - Ρ))\n")
    
    print("Final Answer:")
    print("The time it takes the raindrop to fall is given by the expression:")
    print(f"t = {final_factor} * sqrt(ρ / (ρ - Ρ))")


solve_raindrop_problem()