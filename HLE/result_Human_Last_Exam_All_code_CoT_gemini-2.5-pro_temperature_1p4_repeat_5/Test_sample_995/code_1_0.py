import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

def solve():
    """
    This function implements the step-by-step plan to solve the problem.
    It calculates the minimized probability of Agent A winning and then finds
    the floor of its reciprocal as required.
    """

    # Define the function for the argument of the arccos, which comes from A's optimal strategy.
    # The probability of A winning for r >= F/2 is (1/pi) * arccos(arg).
    def get_arccos_argument(r, F):
        # A's optimal move D is sqrt(2*F*r - F**2). The arccos argument is D/r.
        # The argument for the square root is guaranteed non-negative for r >= F/2.
        sqrt_argument = 2 * F * r - F**2
        return np.sqrt(sqrt_argument) / r

    # Define the function to be integrated.
    # This is P(A wins | r) * p(r), where p(r) = 2r.
    def integrand(r, F):
        arg = get_arccos_argument(r, F)
        
        # Handle potential floating point inaccuracies causing arg > 1.
        if arg > 1.0:
            arg = 1.0
            
        # P(A wins | r) = (1/pi) * arccos(arg)
        # We need to compute Integral[P(A wins | r) * 2r * dr]
        return (1 / np.pi) * np.arccos(arg) * 2 * r

    # Define the total probability of A winning as a function of B's move F.
    def prob_A_wins(F):
        # Agent B's move F must be within the unit circle.
        if F <= 0 or F > 1:
            # Return a large value for invalid F to guide the optimizer.
            return 1.0 

        # The probability calculation is split into two parts based on r vs F/2.
        # Part 1: For r in [0, F/2], A's win probability is 1.
        # The integral is Integral[1 * 2r * dr] from 0 to F/2, which is (F/2)^2.
        prob_part1 = (F**2) / 4

        # Part 2: For r in [F/2, 1], we integrate the numerically calculated probability.
        prob_part2, _ = quad(integrand, F / 2, 1, args=(F,))

        return prob_part1 + prob_part2

    # Use a numerical optimizer to find the value of F that minimizes the win probability.
    # We provide an initial guess (x0) and bounds for F.
    optimization_result = minimize(
        prob_A_wins, 
        x0=0.8, 
        bounds=[(1e-9, 1.0)]
    )

    # The minimized probability of A winning.
    min_prob_A_wins = optimization_result.fun

    # The final value to be computed as per the problem statement.
    final_value = np.floor(1 / min_prob_A_wins)
    
    # Output the components of the final equation as requested.
    print(f"The minimized probability P(A wins) is: {min_prob_A_wins}")
    print(f"The reciprocal 1 / P(A wins) is: {1/min_prob_A_wins}")
    print("The final result is the integer part of the reciprocal.")
    print(f"Final Equation: floor(1 / {min_prob_A_wins:.4f}) = {final_value:.0f}")

solve()