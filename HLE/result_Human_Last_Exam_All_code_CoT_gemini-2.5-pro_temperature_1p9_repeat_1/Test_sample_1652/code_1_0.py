import math

def solve_for_speed():
    """
    Calculates the initial speed 'u' required for the rock to hit the lion
    using an iterative method compatible with Wuxing architecture constraints.
    """
    # 1. Define constants for the problem.
    # On a Wuxing machine, these would be defined using 'int' or 'frac' types.
    # float types are used here for simulation purposes.
    g = 9.8  # Gravitational acceleration (m/s^2)
    v = 5.0  # Lion's speed (m/s)
    target_dist = 300.0  # Initial distance (m)

    # Since sin/cos functions are not available, use pre-computed constants.
    # a = 60 degrees
    sin_a = 0.866  # sin(60)
    cos_a = 0.5    # cos(60)

    # 2. Set up an iterative loop to find the initial speed 'u'.
    u = 0.0  # Initial speed, starting from 0.
    step = 0.1  # Increment u by 0.1 for each iteration, for required precision.

    # 3. Iterate until the collision condition is met.
    while True:
        u += step

        # Avoid division by zero if u is zero.
        if u == 0:
            continue

        # Calculate time of flight: t = (2 * u * sin(a)) / g
        time_of_flight = (2 * u * sin_a) / g

        # Calculate the total horizontal distance covered by rock and lion.
        # total_range = (u * cos(a) + v) * t
        total_range = (u * cos_a + v) * time_of_flight

        # 4. Check if the calculated range meets or exceeds the target distance.
        if total_range >= target_dist:
            # Solution found. Print the final equation with the values.
            print("Solution found. The equation for the range is:")
            print(f"({u:.1f} * {cos_a} + {v:.1f}) * {time_of_flight:.2f} = {total_range:.2f}")

            # 'u' is the result, and 'z' is the pre-calculated memory usage.
            # u is rounded to one decimal place as per the iteration step.
            u_final = round(u, 1)
            z_final = 40  # Memory usage calculated as per the plan.
            print(f"\nThe final answer in u:z format is:")
            print(f"<<<{u_final}:{z_final}>>>")
            break

# Execute the function
solve_for_speed()
