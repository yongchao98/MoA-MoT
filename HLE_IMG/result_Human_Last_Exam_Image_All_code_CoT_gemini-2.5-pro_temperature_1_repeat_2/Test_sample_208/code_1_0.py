def solve_clock_conundrum():
    """
    This function solves the riddle by calculating the meeting time of the two special clock hands.
    """
    # 1. Define speeds in 'minute-markers per minute'.
    # A clock face has 60 minute-markers.
    # The minute hand moves a quarter turn (60 / 4 = 15 markers) per minute.
    speed_minute_hand = 15

    # The hour hand moves in "steps of three", interpreted as 3 markers per minute.
    speed_hour_hand = 3

    # 2. Define the distance the minute hand must gain on the hour hand to meet again.
    # This is one full circle.
    markers_in_a_circle = 60

    # 3. Calculate the relative speed of the minute hand.
    relative_speed = speed_minute_hand - speed_hour_hand

    # 4. Set up and solve the equation for time 't'.
    # The equation is: (speed_minute_hand - speed_hour_hand) * t = markers_in_a_circle
    time_to_meet = markers_in_a_circle / relative_speed

    # 5. Print the explanation and the equation with the numbers.
    print("Solving 'The Clockwork Conundrum':")
    print(f"The minute hand's speed is a quarter turn per minute, which is {speed_minute_hand} markers/min.")
    print(f"The hour hand's speed is 'three steps', which we interpret as {speed_hour_hand} markers/min.")
    print("For the hands to meet, the minute hand must gain one full lap (60 markers) on the hour hand.")
    print("\nThe final equation to find the time in minutes is:")
    print(f"({speed_minute_hand} - {speed_hour_hand}) * time = {markers_in_a_circle}")
    print(f"{relative_speed} * time = {markers_in_a_circle}")
    print(f"time = {markers_in_a_circle} / {relative_speed}")
    print(f"time = {int(time_to_meet)} minutes")

    # 6. Format the answer as a time on a 24-hour clock.
    # The clock starts at noon (12:00).
    start_hour = 12
    meeting_minute = int(time_to_meet)
    final_time = f"{start_hour:02d}:{meeting_minute:02d}"

    print(f"\nSince they start at noon (12:00), they will meet at: {final_time}")

solve_clock_conundrum()