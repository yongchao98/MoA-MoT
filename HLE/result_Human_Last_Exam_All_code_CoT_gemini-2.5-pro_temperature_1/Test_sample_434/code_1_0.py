def find_blowup_condition():
    """
    This function determines the set of initial conditions y(0) for a given x(0) > 1
    that cause the solution of the given system of ODEs to blow up.

    The system is:
    x'(t) = -3*x(t)*y(t)
    y'(t) = -y(t)**2 - x(t) + 1
    """

    # The derivation leads to the analysis of trajectories in the (x, y) plane.
    # The equation for the trajectories is found to be:
    # y**2 = 2*x + 1 + C * x**(2/3)

    # A non-blowing-up solution can converge to the saddle point at (1, 0).
    # The separatrix leading to this point is found by setting x=1, y=0.
    # 0**2 = 2*1 + 1 + C * 1**(2/3) => C = -3.
    # So the separatrix equation is y**2 = 2*x + 1 - 3*x**(2/3).

    # The stable manifold (the trajectory approaching the saddle from x > 1) is
    # y = sqrt(2*x + 1 - 3*x**(2/3)).

    # Solutions blow up if their initial condition (x(0), y(0)) is below this stable manifold.
    # This gives the condition for blow-up.

    # We will now print this condition in a readable format, including all numbers.
    two = 2
    one = 1
    minus_three = -3
    power_two_thirds = "2/3"
    power_one_half = "1/2"

    print("For the system of differential equations:")
    print("x'(t) = -3*x(t)*y(t)")
    print("y'(t) = -y(t)**2 - x(t) + 1")
    print("\nGiven an initial condition (x(0), y(0)) with x(0) > 1, the solution blows up if and only if y(0) satisfies the inequality:")
    
    # Print the final inequality with explicit numbers
    print(f"y(0) < ( {two}*x(0) + {one} - {abs(minus_three)}*x(0)**({power_two_thirds}) )**({power_one_half})")

find_blowup_condition()