import math

def solve_projectile_speed():
    """
    Calculates the initial speed 'u' for the rock to hit the lion and the
    memory usage 'z' based on the Wuxing architecture.
    """
    # Given parameters
    distance_d = 300  # meters
    lion_speed_v = 5  # m/s
    angle_a_deg = 60  # degrees
    gravity_g = 9.8   # m/s^2

    # Convert angle to radians for use in math functions
    angle_a_rad = math.radians(angle_a_deg)

    # The problem can be modeled by the quadratic equation: A*u^2 + B*u + C = 0
    # where u is the initial speed of the rock.
    
    # Calculate the coefficients A, B, and C
    A = math.sin(2 * angle_a_rad)
    B = 2 * lion_speed_v * math.sin(angle_a_rad)
    C = -distance_d * gravity_g

    print("The problem is solved using the quadratic equation: A*u^2 + B*u + C = 0")
    print(f"Calculated coefficients are:")
    print(f"A = sin(2 * {angle_a_deg}°) = {A:.4f}")
    print(f"B = 2 * {lion_speed_v} * sin({angle_a_deg}°) = {B:.4f}")
    print(f"C = -{distance_d} * {gravity_g} = {C:.4f}")
    
    # Solve the quadratic equation for u using the formula: u = [-B + sqrt(B^2 - 4AC)] / 2A
    # We take the positive root because speed cannot be negative.
    discriminant = math.sqrt(B**2 - 4 * A * C)
    initial_rock_speed_u = (-B + discriminant) / (2 * A)
    
    # Round the result to one decimal place
    u_rounded = round(initial_rock_speed_u, 1)

    # Calculate memory usage 'z' in Decimal digits (D) for the Wuxing architecture
    # Variables: distance (int=5D), lion_speed (int=5D), angle (int=5D), 
    # gravity (frac=6D), result_u (frac=6D)
    # Temporary intermediates are excluded as per the instructions.
    mem_distance = 5  # int
    mem_lion_speed = 5 # int
    mem_angle = 5     # int
    mem_gravity = 6   # frac
    mem_result_u = 6  # frac
    z = mem_distance + mem_lion_speed + mem_angle + mem_gravity + mem_result_u

    print("\n---")
    print("Final Answer:")
    # Format the final output as u:z
    print(f"{u_rounded}:{z}")

solve_projectile_speed()
<<<53.5:27>>>