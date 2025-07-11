def solve_hotel_lights():
    """
    Solves the hotel lights puzzle.
    """
    # 0: Red, 1: Green, 2: Blue
    # Use 1-based indexing, so create 101 rooms and ignore index 0.
    num_rooms = 100
    lights = [0] * (num_rooms + 1)

    # Loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        # Guest n toggles every nth room n times.
        for room_number in range(guest_n, num_rooms + 1, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3

        # After each guest, the cat resets any green light to red.
        for room_number in range(1, num_rooms + 1):
            if lights[room_number] == 1:  # 1 is Green
                lights[room_number] = 0   # 0 is Red

    # Find and count the blue lights at the end
    blue_light_rooms = []
    for room_number in range(1, num_rooms + 1):
        if lights[room_number] == 2:  # 2 is Blue
            blue_light_rooms.append(room_number)

    # The problem asks to "output each number in the final equation!".
    # We interpret this as printing the room numbers that end up blue.
    if blue_light_rooms:
        print("The following room numbers have blue lights:")
        # The equation can be seen as the sum of these room numbers
        equation = " + ".join(map(str, blue_light_rooms))
        print(f"{equation} = {sum(blue_light_rooms)}")
    else:
        print("No rooms have blue lights.")

    print("\nThe total number of blue lights is:")
    print(len(blue_light_rooms))

solve_hotel_lights()
<<<33>>>