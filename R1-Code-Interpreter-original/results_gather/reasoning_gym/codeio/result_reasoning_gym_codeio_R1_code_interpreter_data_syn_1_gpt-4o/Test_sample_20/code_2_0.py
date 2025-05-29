def find_feasible_input():
    # Given output
    target_throttle_needed = -2.801755499249228
    target_should_start_burn = True

    # Define broader ranges for inputs with smaller increments
    max_thrust_range = [i for i in range(500, 2000, 100)]  # Newtons
    current_mass_range = [i for i in range(400, 1000, 50)]  # kg
    speed_diff_range = [i for i in range(30, 70, 5)]  # m/s
    time_to_impact_range = [i for i in range(5, 15, 1)]  # seconds
    drag_range = [i for i in range(100, 300, 10)]  # Newtons
    gravitational_accel = 9.81  # m/s^2
    current_speed_range = [i for i in range(50, 150, 10)]  # m/s

    # Iterate over all combinations of input ranges
    for max_thrust in max_thrust_range:
        for current_mass in current_mass_range:
            for speed_diff in speed_diff_range:
                for time_to_impact in time_to_impact_range:
                    for drag in drag_range:
                        for current_speed in current_speed_range:
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
    return None

# Execute the function to find feasible input
print(find_feasible_input())