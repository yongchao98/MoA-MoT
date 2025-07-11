import fractions

def solve_particle_problem():
    """
    Solves the random walk problem by analyzing the steady-state of relative gaps.
    """
    # Step 1: Solve the system of linear equations for gap probabilities.
    # The system is:
    #   2*p1 - p2 = 1/3
    #   -p1 + 2*p2 = 1
    # We solve this system for p1 = P(Y1 > 1) and p2 = P(Y2 > 1).
    # From the second eq: p1 = 2*p2 - 1.
    # Substitute into the first: 2*(2*p2 - 1) - p2 = 1/3  => 3*p2 = 7/3 => p2 = 7/9
    # Then p1 = 2*(7/9) - 1 = 14/9 - 9/9 = 5/9
    p1 = fractions.Fraction(5, 9)
    p2 = fractions.Fraction(7, 9)

    print("--- Step 1: Solving for Gap Probabilities ---")
    print("Let p1 = P(Y1 > 1) and p2 = P(Y2 > 1). The steady-state balance equations are:")
    print("  2*p1 - 1*p2 = 1/3")
    print(" -1*p1 + 2*p2 = 1")
    print(f"The solution is p1 = {p1} and p2 = {p2}\n")

    # Step 2: Calculate average gap sizes.
    # The gaps Y1 and Y2 follow geometric distributions with success probabilities (1-rho).
    # The parameter rho is equal to the probability of the gap being > 1.
    # rho1 = p1, rho2 = p2.
    # The mean of a shifted geometric distribution (starting at 1) is 1/(1-rho).
    rho1 = p1
    rho2 = p2
    mean_y1 = 1 / (1 - rho1)
    mean_y2 = 1 / (1 - rho2)

    print("--- Step 2: Calculating Average Gap Distances ---")
    print("The gaps Y1 and Y2 follow independent geometric distributions.")
    print(f"The average size of gap Y1 is calculated as:")
    print(f"  E[Y1] = 1 / (1 - {rho1}) = {mean_y1}")
    print(f"The average size of gap Y2 is calculated as:")
    print(f"  E[Y2] = 1 / (1 - {rho2}) = {mean_y2}\n")

    # Step 3: Calculate the total average distance.
    total_distance = mean_y1 + mean_y2
    print("--- Step 3: Calculating Total Average Distance ---")
    print("The total average distance is the sum of the average gap sizes:")
    print(f"  E[Distance] = E[Y1] + E[Y2] = {mean_y1} + {mean_y2} = {total_distance}\n")

    # Step 4: Calculate the asymptotic speed.
    # The speed is the average velocity of any particle. For the leftmost particle:
    # v = (jump rate right) * (+1) + (jump rate left) * (-1)
    # The jump right is conditioned on Y1 > 1. The jump left is unconditional.
    rate_right_p1 = p1
    rate_left_p1 = fractions.Fraction(1, 3)
    speed = rate_right_p1 - rate_left_p1

    print("--- Step 4: Calculating Asymptotic Speed ---")
    print("The asymptotic speed is the long-term average velocity of the system.")
    print("Calculated for the leftmost particle:")
    print(f"  v = P(Y1 > 1) * (+1) + (1/3) * (-1)")
    print(f"  v = {rate_right_p1} - {rate_left_p1} = {speed}\n")

    # Step 5: Final Result.
    print("--- Final Answer ---")
    print("The average distance and asymptotic speed are:")
    print(f"({total_distance}, {speed})")

solve_particle_problem()