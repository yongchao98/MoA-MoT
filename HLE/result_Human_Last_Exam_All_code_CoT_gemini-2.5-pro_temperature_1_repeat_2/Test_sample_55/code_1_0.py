def solve_hotel_lights():
    """
    Solves the hotel lights riddle by simulating the actions of guests and a cat.
    """
    # 0: Red, 1: Green, 2: Blue
    # Initialize 100 rooms (indices 1-100), all lights are red (0).
    lights = [0] * 101

    # Loop through guests from 1 to 100
    for guest_n in range(1, 101):
        # Guest n toggles the light in every nth room n times.
        for room_i in range(guest_n, 101, guest_n):
            lights[room_i] = (lights[room_i] + guest_n) % 3

        # After each guest leaves, a cat resets any green light to red.
        for room_i in range(1, 101):
            if lights[room_i] == 1:  # If Green
                lights[room_i] = 0  # Reset to Red

    # Find all rooms with blue lights at the end.
    blue_light_rooms = []
    for room_i in range(1, 101):
        if lights[room_i] == 2:  # If Blue
            blue_light_rooms.append(room_i)

    # As requested, output the numbers that contribute to the final answer.
    # First, we print the room numbers that are blue.
    print("The following rooms have a blue light:")
    print(blue_light_rooms)

    # Then, we print the final count.
    count = len(blue_light_rooms)
    print(f"\nTotal number of blue lights: {count}")

solve_hotel_lights()
<<<33>>>