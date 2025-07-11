def solve_clock_conundrum():
    """
    Solves the Clockwork Conundrum riddle.
    """
    # 1. Define the clock and hand speeds based on the riddle.
    # The clock has 60 markers (like minutes).
    TOTAL_MARKERS = 60
    # Hour hand speed: 3 markers per tick.
    hour_hand_speed = 3
    # Minute hand speed: A quarter turn (1/4 of 60) is 15 markers per tick.
    minute_hand_speed = 15

    # 2. Calculate the relative speed.
    relative_speed = minute_hand_speed - hour_hand_speed

    # 3. Calculate how many ticks until they meet.
    # This happens when the minute hand has lapped the hour hand by one full circle.
    # Time = Total Distance / Relative Speed
    ticks_to_meet = TOTAL_MARKERS // relative_speed

    # 4. Find the exact marker on the clock where they meet.
    # We can use either hand's speed for this.
    meeting_position = hour_hand_speed * ticks_to_meet

    # 5. Convert the meeting position to a standard time (HH:MM).
    # The hour is the clock number (position / 5).
    # The minute is simply the position marker.
    hour = meeting_position // 5
    minute = meeting_position
    
    # As requested, here are the numbers used in the final equation.
    print("Equation to find meeting position:")
    print(f"{hour_hand_speed} * ({TOTAL_MARKERS} / ({minute_hand_speed} - {hour_hand_speed})) = {meeting_position}")
    
    print("\nEquation for the final time:")
    print(f"Hour = {meeting_position} / 5 = {hour}")
    print(f"Minute = {meeting_position}")

    print("\nFinal time:")
    # Print the final time in HH:MM format, respecting airline convention style.
    print(f"{hour}:{minute:02d}")

solve_clock_conundrum()