def solve_hotel_lights_riddle():
    """
    Solves the hotel lights riddle by simulating the actions of 100 guests and a cat.
    """
    # Define the states for the lights: 0 for Red, 1 for Green, 2 for Blue.
    # We use a list of size 101 so we can use 1-based indexing for the 100 rooms.
    num_rooms = 100
    lights = [0] * (num_rooms + 1) # All lights start as Red (0).

    # Loop through each guest from 1 to 100.
    for guest_n in range(1, num_rooms + 1):
        # Guest 'n' toggles the light 'n' times in every n-th room.
        for room_k in range(guest_n, num_rooms + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3

        # After each guest finishes, the cat resets all Green lights to Red.
        for room_k in range(1, num_rooms + 1):
            if lights[room_k] == 1:  # If the light is Green...
                lights[room_k] = 0   # ...the cat resets it to Red.

    # After all guests have left, find which rooms have blue lights.
    blue_rooms = []
    for room_k in range(1, num_rooms + 1):
        if lights[room_k] == 2:  # If the light is Blue...
            blue_rooms.append(room_k)

    # Print the final result.
    num_blue_lights = len(blue_rooms)
    print(f"The following {num_blue_lights} rooms have blue lights:")
    print(blue_rooms)
    
    # As requested, output each number in the final equation for the total count.
    if num_blue_lights > 0:
        equation_parts = ['1'] * num_blue_lights
        equation_str = " + ".join(equation_parts)
        print("\nThe final equation to find the total number is:")
        print(f"{equation_str} = {num_blue_lights}")
    else:
        print("\nNo lights are blue, so the total count is 0.")

# Run the simulation and print the results.
solve_hotel_lights_riddle()
<<<33>>>