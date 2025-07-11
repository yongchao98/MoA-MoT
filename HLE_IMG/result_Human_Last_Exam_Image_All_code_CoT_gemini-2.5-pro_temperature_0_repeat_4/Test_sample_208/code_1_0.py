def solve_clock_conundrum():
    """
    Solves the riddle by calculating the meeting time of the two strange clock hands.
    """
    # 1. Define the parameters from the riddle.
    # A standard clock face has 12 hour-positions.
    clock_face_positions = 12

    # The hour hand's speed is 3 positions per tick.
    hour_hand_speed = 3

    # The minute hand's relative speed ("gains") is a quarter turn per tick.
    # A quarter turn on a 12-position clock is 12 / 4 = 3 positions.
    relative_speed = clock_face_positions / 4

    # 2. Calculate the number of ticks until they meet.
    # They meet when the relative distance covered is one full circle.
    # Time = Distance / Speed
    ticks_to_meet = clock_face_positions / relative_speed

    # 3. Determine the final time.
    # The clock starts at noon (12:00). We assume a "tick" is one hour.
    # In 24-hour format, noon is 12.
    start_hour = 12
    meeting_hour = start_hour + ticks_to_meet

    # As requested, output each number in the final equations.
    print("First, we calculate the number of ticks (hours) until the hands meet again.")
    print("The equation is: Total Positions / Relative Speed = Ticks to Meet")
    print(f"{int(clock_face_positions)} / {int(relative_speed)} = {int(ticks_to_meet)}")
    
    print("\nNext, we calculate the time of day in 24-hour format.")
    print("The equation is: Start Hour + Hours Elapsed = Meeting Hour")
    print(f"{int(start_hour)} + {int(ticks_to_meet)} = {int(meeting_hour)}")

    # 4. Format the final answer.
    # The meeting happens at an exact hour, so the minutes are 00.
    final_time_str = f"{int(meeting_hour)}:00"
    
    print(f"\nTherefore, the hands will meet again at: {final_time_str}")

solve_clock_conundrum()