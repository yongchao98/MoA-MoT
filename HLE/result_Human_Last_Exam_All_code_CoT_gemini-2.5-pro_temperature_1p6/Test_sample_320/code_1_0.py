from fractions import Fraction

def solve_particle_problem():
    """
    Solves for the average distance and asymptotic speed of a 3-particle system.

    The problem is modeled by analyzing the gaps between particles, Y1 and Y2.
    The stationary distribution of these gaps is found to be a product of two
    geometric distributions with parameters alpha and beta.
    pi(i,j) ~ alpha**(i-1) * beta**(j-1)

    The parameters alpha and beta are solved from the balance equations at the
    boundaries of the state space.
    - At state (1,1): alpha + beta = 4/3
    - At state (i,1), i>1: 10*alpha/3 = 1/3 + alpha**2 + beta*(1+alpha)

    Substituting beta = 4/3 - alpha into the second equation:
    10*alpha/3 = 1/3 + alpha**2 + (4/3 - alpha)*(1+alpha)
    10*alpha/3 = 1/3 + alpha**2 + 4/3 + 4*alpha/3 - alpha - alpha**2
    10*alpha/3 = 5/3 + alpha/3
    9*alpha/3 = 5/3
    3*alpha = 5/3
    alpha = 5/9
    """

    # Rates for the leftmost particle
    p_L = Fraction(1, 3)
    p_R = 1

    # Solve for alpha and beta
    alpha = Fraction(5, 9)
    beta = Fraction(4, 3) - alpha

    print(f"Solved parameters for the stationary distribution:")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}\n")

    # The distributions for Y1 and Y2 are geometric.
    # P(Y1=i) = (1-alpha) * alpha**(i-1)
    # P(Y2=j) = (1-beta) * beta**(j-1)
    # The mean of such a distribution is 1/(1-r).

    E_Y1 = 1 / (1 - alpha)
    E_Y2 = 1 / (1 - beta)

    print("Calculating the average gaps:")
    print(f"E[Y1] = 1 / (1 - {alpha}) = {E_Y1}")
    print(f"E[Y2] = 1 / (1 - {beta}) = {E_Y2}\n")

    # Calculate the average distance between the leftmost and rightmost particles
    avg_distance = E_Y1 + E_Y2
    
    print("Calculating the average distance between the leftmost and rightmost particles:")
    print(f"D = E[Y1] + E[Y2] = {E_Y1} + {E_Y2} = {avg_distance}\n")

    # Calculate the asymptotic speed of the leftmost particle
    # v1 = (+1)*p_R*P(Y1>1) + (-1)*p_L
    # P(Y1>1) = 1 - P(Y1=1) = 1 - (1-alpha) = alpha
    speed = p_R * alpha - p_L

    print("Calculating the asymptotic speed of the leftmost particle:")
    print(f"v1 = (rate right) * P(Y1>1) - (rate left)")
    print(f"v1 = {p_R} * {alpha} - {p_L} = {speed}\n")

    # Final result
    result_tuple = (avg_distance, speed)
    print(f"The final result (distance, speed) is: {result_tuple}")

if __name__ == '__main__':
    solve_particle_problem()