def solve_light_puzzle():
    """
    Solves the hotel light puzzle and prints the result.
    Colors are represented as: 0=Red, 1=Green, 2=Blue.
    """
    # 1. Initialize 100 rooms (plus a placeholder at index 0).
    # All lights are initially Red (0).
    lights = [0] * 101

    # 2. Loop through each guest from 1 to 100.
    for guest_n in range(1, 101):
        # 3. Guest n toggles every nth room n times.
        for room_number in range(guest_n, 101, guest_n):
            # Toggling n times means adding n to the current color state, modulo 3.
            lights[room_number] = (lights[room_number] + guest_n) % 3

        # 4. After each guest, the cat resets any green light to red.
        for i in range(1, 101):
            if lights[i] == 1:  # If the light is Green
                lights[i] = 0   # Reset to Red

    # 5. Find the final count and the room numbers of blue lights.
    blue_lights_rooms = []
    for i in range(1, 101):
        if lights[i] == 2:  # If the light is Blue
            blue_lights_rooms.append(i)

    # 6. Output the results.
    print(f"The number of lights that will be blue is: {len(blue_lights_rooms)}")
    print("The equation for the final blue lights is the list of room numbers that are blue:")
    
    # As requested, output each number in the final equation (the list of rooms)
    # The 'equation' is simply the sum of the blue lights, or the list of them.
    # We will print the list in a format resembling an equation.
    if not blue_lights_rooms:
        print("0 = 0")
    else:
        # Construct the string like "1 + 1 + ... = total"
        count_str = " + ".join(["1"] * len(blue_lights_rooms))
        print(f"The rooms that are blue are: {', '.join(map(str, blue_lights_rooms))}")
        print(f"Representing each blue light as '1', the equation is: {count_str} = {len(blue_lights_rooms)}")


solve_light_puzzle()
<<<20>>>