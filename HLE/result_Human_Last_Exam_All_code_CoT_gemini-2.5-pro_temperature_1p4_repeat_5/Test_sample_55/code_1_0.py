def solve_light_riddle():
    """
    Solves the hotel light riddle.

    - 100 rooms, lights cycle red -> green -> blue (0 -> 1 -> 2).
    - All lights start red (0).
    - 100 guests arrive. Guest 'n' toggles every nth room 'n' times.
    - Toggling 'n' times means new_state = (current_state + n) % 3.
    - After each guest, a cat resets any green (1) light to red (0).
    - The function determines the number of blue lights at the end.
    """

    # 0: red, 1: green, 2: blue
    # Use list of size 101 to have indices 1-100 for rooms 1-100.
    lights = [0] * 101

    # Loop for each guest from 1 to 100
    for n in range(1, 101):
        # Guest 'n' visits every nth room
        for room_number in range(n, 101, n):
            # Toggle the light 'n' times
            lights[room_number] = (lights[room_number] + n) % 3

        # After the guest leaves, the cat resets green lights to red
        for i in range(1, 101):
            if lights[i] == 1:  # If light is green
                lights[i] = 0   # Reset to red

    # After all guests have left, find the rooms with blue lights
    blue_light_rooms = []
    for i in range(1, 101):
        if lights[i] == 2:  # If light is blue
            blue_light_rooms.append(i)

    # Print the "equation" as requested, showing each room number
    if not blue_light_rooms:
        print("There are no blue lights at the end.")
        print("Total number of blue lights: 0")
    else:
        # Constructing the equation "1 + 1 + ... = count" visually
        # by listing the room numbers that contribute to the count.
        print("The following room numbers have blue lights:")
        equation_parts = [str(room) for room in blue_light_rooms]
        print(' '.join(equation_parts))
        print(f"Total number of blue lights: {len(blue_light_rooms)}")

solve_light_riddle()
<<<20>>>