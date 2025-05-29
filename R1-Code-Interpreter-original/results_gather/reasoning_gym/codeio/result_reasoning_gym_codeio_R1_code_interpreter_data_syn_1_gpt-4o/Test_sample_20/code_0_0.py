def find_feasible_input():
    # Given output
    target_throttle_needed = -2.801755499249228
    target_should_start_burn = True

    # Initial guesses for inputs
    max_thrust = 1000.0  # Newtons
    current_mass = 500.0  # kg
    speed_diff = 50.0  # m/s
    time_to_impact = 10.0  # seconds
    drag = 200.0  # Newtons
    gravitational_accel = 9.81  # m/s^2
    current_speed = 100.0  # m/s

    # Calculate needed_force for throttle_needed
    desired_accel = speed_diff / time_to_impact
    accel_force = current_mass * desired_accel
    current_force = (current_mass * gravitational_accel * -1) + drag
    needed_force = accel_force + current_force
    throttle_needed = needed_force / max_thrust

    # Calculate needed_force for should_start_burn
    inertial_force = current_mass * (current_speed / time_to_impact)
    needed_force_burn = inertial_force + current_force - drag
    should_start_burn = max_thrust <= needed_force_burn

    # Check if the calculated values match the target output
    if abs(throttle_needed - target_throttle_needed) < 0.01 and should_start_burn == target_should_start_burn:
        return {
            "max_thrust": max_thrust,
            "current_mass": current_mass,
            "speed_diff": speed_diff,
            "time_to_impact": time_to_impact,
            "drag": drag,
            "gravitational_accel": gravitational_accel,
            "current_speed": current_speed
        }
    else:
        return None

# Execute the function to find feasible input
print(find_feasible_input())