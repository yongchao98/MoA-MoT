import math

def solve_clock_conundrum():
    """
    Solves the Clockwork Conundrum riddle by calculating the meeting time of the two hands.
    """
    
    # 1. Define constants based on a standard clock.
    DEGREES_IN_CIRCLE = 360

    # 2. Analyze the hour hand's movement.
    # The riddle states "The hour hand moves in steps of three".
    # This implies its speed is 3 times the normal hour hand speed.
    # A normal hour hand moves 360 degrees in 12 hours.
    # Its speed is (360 / 12) = 30 degrees per hour.
    # The conundrum's hour hand speed is 3 * 30 = 90 degrees per hour.
    hour_hand_speed_deg_per_hr = (DEGREES_IN_CIRCLE / 12) * 3

    # 3. Analyze the minute hand's movement.
    # The riddle states "The minute hand gains, with every tick, a quarter turn".
    # A "quarter turn" is 360 / 4 = 90 degrees.
    # A "tick" in the context of a clock is best interpreted as a minute.
    # The keyword "gains" implies this is its speed RELATIVE to the hour hand.
    relative_speed_deg_per_min = DEGREES_IN_CIRCLE / 4

    # 4. Formulate the equation.
    # The hands start together. They meet when the minute hand has "lapped" the
    # hour hand, which means it has gained a full 360 degrees on it.
    # Time = Total Angle to Gain / Relative Speed
    # We need to make sure the units are consistent (degrees per minute).
    
    # Let's verify our units and numbers for the equation.
    angle_to_gain = DEGREES_IN_CIRCLE
    # The relative speed is already defined in degrees per minute.
    relative_speed = relative_speed_deg_per_min

    # 5. Solve for the time in minutes.
    time_to_meet_min = angle_to_gain / relative_speed
    
    # Print out the thinking process and the final equation.
    print("The plan is to find the time it takes for the minute hand to gain 360 degrees on the hour hand.")
    print(f"The relative speed (gain) of the minute hand is a 'quarter turn' per minute, which is {int(relative_speed)} degrees per minute.")
    print("The total angle to gain for the hands to meet again is a full circle.")
    print(f"\nThe final equation is:")
    print(f"Time (minutes) = {int(angle_to_gain)} degrees / {int(relative_speed)} degrees per minute")
    
    result_minutes = int(time_to_meet_min)
    print(f"Result = {result_minutes} minutes")

    # 6. Format the final answer.
    # The start time is noon (12:00).
    # The meeting time is 'result_minutes' past noon.
    start_hour = 12
    start_minute = 0
    
    final_minute = start_minute + result_minutes
    final_hour = start_hour 

    # Format with leading zero for minutes if needed.
    final_time_str = f"{final_hour:02d}:{final_minute:02d}"
    
    print("\nThe hands will meet at the following time:")
    print(final_time_str)

solve_clock_conundrum()