def find_feasible_input():
    # Further adjusted guesses for the inputs
    max_thrust = 3000.0  # Newtons
    current_mass = 900.0  # kg
    speed_diff = -25.0  # m/s
    time_to_impact = 18.0  # seconds
    drag = 250.0  # Newtons
    gravitational_accel = 9.81  # m/s^2
    current_speed = 200.0  # m/s

    # Calculate current_force
    current_force = (current_mass * gravitational_accel * -1) + drag

    # Calculate desired_accel
    desired_accel = speed_diff / time_to_impact

    # Calculate accel_force
    accel_force = current_mass * desired_accel

    # Calculate needed_force for throttle
    needed_force_throttle = accel_force + current_force

    # Calculate throttle_needed
    throttle_needed = needed_force_throttle / max_thrust

    # Calculate inertial_force
    inertial_force = current_mass * (current_speed / time_to_impact)

    # Calculate needed_force for burn
    needed_force_burn = inertial_force + current_force - drag

    # Determine should_start_burn
    should_start_burn = max_thrust <= needed_force_burn

    return {
        "max_thrust": max_thrust,
        "current_mass": current_mass,
        "speed_diff": speed_diff,
        "time_to_impact": time_to_impact,
        "drag": drag,
        "gravitational_accel": gravitational_accel,
        "current_speed": current_speed,
        "throttle_needed": throttle_needed,
        "should_start_burn": should_start_burn
    }

print(find_feasible_input())