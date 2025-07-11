import numpy as np
from fractions import Fraction

def solve_particle_walk():
    """
    Calculates the average distance and asymptotic speed for the described
    three-particle system.
    """
    # Step 1: Explain the methodology and setup
    print("This problem describes a three-particle continuous-time random walk with exclusion.")
    print("The system is a variant of the Asymmetric Simple Exclusion Process (ASEP).")
    print("We analyze the system by looking at the gaps between particles: y1 = x2 - x1 and y2 = x3 - x2.")
    print("The stationary state distribution for these gaps is assumed to be a product of geometric distributions:")
    print("P(y1, y2) = C * alpha**(y1-1) * beta**(y2-1)")
    print("\nNext, we express the average velocity 'v' for each particle.")
    print("The velocity depends on the jump rates and the probability of the target site being empty.")
    print("Let p_i and q_i be the right and left jump rates for particle i.")
    print("p1=1, q1=1/3; p2=1, q2=1; p3=1, q3=1.")
    print("The probability that a particle can jump is related to the gaps:")
    print("P(site x_i+1 is empty) = P(y_i > 1), P(site x_i-1 is empty) = P(y_{i-1} > 1).")
    print("For the geometric distribution, P(y > 1) = parameter.")
    print("v1 = p1*P(y1 > 1) - q1 = 1*alpha - 1/3")
    print("v2 = p2*P(y2 > 1) - q2*P(y1 > 1) = 1*beta - 1*alpha")
    print("v3 = p3 - q3*P(y2 > 1) = 1 - 1*beta")
    
    # Step 2: Set up and solve the system of equations
    print("\nIn steady state, all particles must have the same speed: v1 = v2 = v3.")
    print("This gives a system of linear equations for alpha and beta:")
    print("1) alpha - 1/3 = beta - alpha  =>  2*alpha - beta = 1/3")
    print("2) beta - alpha = 1 - beta    =>  alpha - 2*beta = -1")

    # We solve this system manually for precision with fractions.
    # From (2), alpha = 2*beta - 1.
    # Substitute into (1): 2*(2*beta - 1) - beta = 1/3
    # 4*beta - 2 - beta = 1/3  =>  3*beta = 7/3  => beta = 7/9
    # alpha = 2*(7/9) - 1 = 14/9 - 9/9 = 5/9
    alpha = Fraction(5, 9)
    beta = Fraction(7, 9)
    
    print("\nSolving this system yields:")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}")

    # Step 3: Calculate the average distance
    print("\nThe average distance between the leftmost and rightmost particles is E[x3 - x1] = E[y1] + E[y2].")
    print("For a geometric distribution, the mean is E[y] = 1 / (1 - parameter).")
    
    E_y1 = 1 / (1 - alpha)
    E_y2 = 1 / (1 - beta)
    E_D = E_y1 + E_y2

    print(f"E[y1] = 1 / (1 - {alpha}) = 1 / ({1-alpha}) = {E_y1}")
    print(f"E[y2] = 1 / (1 - {beta}) = 1 / ({1-beta}) = {E_y2}")
    print(f"Total Average Distance = {E_y1} + {E_y2} = {E_D}")

    # Step 4: Calculate the asymptotic speed
    print("\nThe asymptotic speed 'v' is calculated from any of the velocity equations.")
    print("Using v = 1 - beta:")
    
    v = 1 - beta
    
    print(f"v = 1 - {beta} = {v}")

    # Step 5: Final Answer
    print("\nThe final result (distance, speed) is:")
    print(f"({E_D}, {v})")

solve_particle_walk()