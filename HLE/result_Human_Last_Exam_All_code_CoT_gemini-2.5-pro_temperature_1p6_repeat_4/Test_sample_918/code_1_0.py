import math

def solve_moth_displacement():
    """
    Calculates the moth's displacement when the easternmost LED blinks.
    """
    # Define problem parameters
    moth_airspeed_mpm = 5.675  # Moth's speed in m/min
    wind_speed_mpm = 5.0      # Wind speed in m/min (eastward, against the moth)
    tunnel_length_m = 2.0     # Tunnel length in meters
    num_leds = 80             # Total number of LEDs
    led_blink_delay_s = 0.3   # Time delay between consecutive LED blinks in seconds
    
    # The moth starts at the eastern end and flies west.
    # The blinking sequence starts when the moth is halfway through the tunnel.
    # Displacement at the start of the blinking sequence:
    start_of_blinking_pos_m = tunnel_length_m / 2.0
    
    # Calculate the moth's ground speed in m/min
    # This is the effective speed of the moth relative to the ground.
    ground_speed_mpm = moth_airspeed_mpm - wind_speed_mpm
    
    # Calculate the total time it takes for the LED signal to propagate from the
    # westernmost LED to the easternmost LED. There are (num_leds - 1) intervals.
    led_wave_time_s = (num_leds - 1) * led_blink_delay_s
    
    # Convert the propagation time from seconds to minutes to match the speed units.
    led_wave_time_min = led_wave_time_s / 60.0
    
    # Calculate the additional distance the moth travels during this propagation time.
    additional_distance_m = ground_speed_mpm * led_wave_time_min
    
    # The final displacement is the position where the blinking started (halfway point)
    # plus the additional distance traveled.
    final_displacement_m = start_of_blinking_pos_m + additional_distance_m
    
    print("The moth's final displacement is its position when the blinking starts (halfway point) plus the additional distance it travels while the LED signal propagates from west to east.")
    print("\nThe equation for the final displacement is:")
    print("Final Displacement = Halfway Distance + (Moth Ground Speed * LED Propagation Time)")
    
    print("\nBreaking down the calculation:")
    # We display the final equation using the original numbers from the problem statement.
    # Halfway Distance is 1.0 m.
    # Moth Ground Speed = (moth_airspeed_mpm - wind_speed_mpm)
    # LED Propagation Time in minutes = ((num_leds - 1) * led_blink_delay_s) / 60
    print(f"Final Displacement = {start_of_blinking_pos_m}m + ({moth_airspeed_mpm}m/min - {wind_speed_mpm}m/min) * (({num_leds} - 1) * {led_blink_delay_s}s / 60s/min)")
    print(f"Final Displacement = {start_of_blinking_pos_m} + ({ground_speed_mpm}) * ({led_wave_time_min})")
    print(f"Final Displacement = {final_displacement_m} m")
    
    print("\nThis displacement of ~1.27m corresponds to answer choice D.")

solve_moth_displacement()
<<<D>>>