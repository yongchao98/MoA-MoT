def solve_clock_conundrum():
    """
    This function solves the clock riddle by simulating the movement of the hands.

    It interprets the riddle's rules as follows:
    - The minute hand moves 90 degrees every minute.
    - The hour hand jumps 90 degrees every 4 minutes.

    The function finds the first time 't' (in minutes past noon) where the
    angular positions of both hands are identical.
    """

    time_in_minutes = 0
    while True:
        # We start checking from 1 minute past noon.
        time_in_minutes += 1

        # Calculate the minute hand's position in degrees.
        # It moves a quarter turn (90 degrees) per minute.
        minute_hand_pos = (time_in_minutes * 90) % 360

        # Calculate the hour hand's position in degrees.
        # It moves in a step of 3 hours (90 degrees) every 4 minutes.
        # The number of jumps is the number of 4-minute intervals passed.
        num_jumps = time_in_minutes // 4
        hour_hand_pos = (num_jumps * 90) % 360

        # Check if the hands are at the same position.
        if minute_hand_pos == hour_hand_pos:
            meeting_minutes = time_in_minutes
            break

    # The clock starts at noon (12:00). The hands meet after 'meeting_minutes'.
    # As per the riddle's persona (Qatar Airways convention), we use 24-hour format.
    hour = 12
    minute = meeting_minutes
    
    # The prompt requests an answer with no letters, only numbers and colons.
    # It also says to "output each number in the final equation!".
    # We interpret this as printing the final time, which contains the numbers for the hour and minute.
    # The final equation is Time = Hour : Minute
    print(f"{hour}:{minute:02d}")

solve_clock_conundrum()