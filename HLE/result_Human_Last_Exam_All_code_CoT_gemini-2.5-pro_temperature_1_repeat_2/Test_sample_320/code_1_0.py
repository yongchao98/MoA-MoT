from fractions import Fraction

def solve_particle_problem():
    """
    Solves for the average distance and asymptotic speed of a three-particle system.
    """
    
    # --- Introduction ---
    print("This script calculates the average distance and asymptotic speed for a system of three interacting particles on a line.")
    print("The problem is solved by analyzing the stationary state of the gaps between particles.\n")
    
    # --- Step 1: Solve for stationary probabilities ---
    print("Step 1: Determine the probabilities of minimal gaps.")
    print("Let p1 = P(gap1=1) and p2 = P(gap2=1).")
    print("From the stationary state condition, we derive the following linear system:")
    print("2*p1 - p2 = 2/3")
    print("-p1 + 2*p2 = 0\n")

    # Solving the system:
    # From -p1 + 2*p2 = 0  => p1 = 2*p2
    # Substitute into the first equation: 2*(2*p2) - p2 = 2/3 => 4*p2 - p2 = 2/3 => 3*p2 = 2/3
    # This gives p2 = 2/9.
    # And p1 = 2 * (2/9) = 4/9.
    
    p1 = Fraction(4, 9)
    p2 = Fraction(2, 9)
    
    print(f"Solving this system yields:")
    print(f"p1 = P(Y1=1) = {p1}")
    print(f"p2 = P(Y2=1) = {p2}\n")

    # --- Step 2: Calculate average distance ---
    print("Step 2: Calculate the average distance between the leftmost and rightmost particles.")
    # Assuming a geometric distribution for the gaps, the mean gap size E[Y] is 1/P(Y=1).
    E_Y1 = 1 / p1
    E_Y2 = 1 / p2
    total_distance = E_Y1 + E_Y2
    
    print("The total average distance is the sum of the average gaps, E[Y1] + E[Y2].")
    print(f"The final equation for the distance is: {E_Y1} + {E_Y2} = {total_distance}\n")

    # --- Step 3: Calculate asymptotic speed ---
    print("Step 3: Calculate the asymptotic speed of the leftmost particle.")
    # The speed v1 is the net current: (rate right) * P(jump possible) - (rate left).
    # Jump right is possible if the gap Y1 > 1. P(Y1 > 1) = 1 - p1.
    prob_Y1_gt_1 = 1 - p1
    rate_left = Fraction(1, 3)
    speed = prob_Y1_gt_1 * 1 - rate_left
    
    print("The speed is v1 = P(Y1 > 1) * (rate right) - (rate left).")
    print(f"The final equation for the speed is: ({prob_Y1_gt_1})*1 - {rate_left} = {speed}\n")

    # --- Step 4: Final Answer ---
    print("The calculated quantities are:")
    print(f"Average distance = {total_distance}")
    print(f"Asymptotic speed = {speed}")
    print("\nFinal answer in the format (distance, speed):")
    print(f"({total_distance}, {speed})")

if __name__ == '__main__':
    solve_particle_problem()
