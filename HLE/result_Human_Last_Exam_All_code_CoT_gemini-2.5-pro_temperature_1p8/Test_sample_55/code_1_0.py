def solve_light_puzzle():
    """
    Solves the hotel light puzzle by simulating the actions of 100 guests and a cat.
    """
    # Step 1: Initialize 100 rooms. We use a list of size 101 for 1-based indexing.
    # States: 0 = Red, 1 = Green, 2 = Blue. All lights start as Red.
    lights = [0] * 101

    # Step 2: Loop through each of the 100 guests.
    for guest_n in range(1, 101):
        
        # Step 2a: Guest 'n' visits every nth room and toggles the light 'n' times.
        for room_number in range(guest_n, 101, guest_n):
            current_state = lights[room_number]
            # Toggling 'n' times is equivalent to adding 'n' to the state, modulo 3.
            lights[room_number] = (current_state + guest_n) % 3
            
        # Step 2b: After each guest leaves, a cat resets any Green light to Red.
        for i in range(1, 101):
            if lights[i] == 1:  # If light is Green
                lights[i] = 0   # Reset to Red

    # Step 3: After all guests have left, find all rooms with a Blue light.
    blue_light_rooms = []
    for room_number in range(1, 101):
        if lights[room_number] == 2:  # If light is Blue
            blue_light_rooms.append(room_number)

    # Step 4: Prepare the final output string and print the result.
    count = len(blue_light_rooms)
    if count > 0:
        equation_str = " + ".join(map(str, blue_light_rooms))
        print(f"The rooms with blue lights are numbers: {equation_str}")
        print(f"Adding these room numbers is not required, the count is the answer.")
        print(f"The equation representing the rooms that are blue is effectively the list itself, leading to a total count.")
        # Final output of the sum requested by the prompt.
        print(f"{equation_str} leads to a total of {count} blue lights.")
    else:
        print("There are no blue lights at the end.")
        
    print(f"\nThe final count of blue lights is: {count}")

solve_light_puzzle()
<<<32>>>