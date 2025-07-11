def solve_light_puzzle():
    """
    Solves the hotel light puzzle and prints the result.
    """
    # Define states for clarity
    RED = 0
    GREEN = 1
    BLUE = 2

    NUM_ROOMS = 100

    # Initialize all 100 rooms with red lights.
    # Use a list of size 101 for 1-based indexing (ignore index 0).
    lights = [RED] * (NUM_ROOMS + 1)

    # Loop through each of the 100 guests
    for guest_n in range(1, NUM_ROOMS + 1):
        # Guest 'n' toggles every 'n'th room 'n' times.
        for room_num in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_num] = (lights[room_num] + guest_n) % 3

        # After each guest leaves, the cat resets any green lights to red.
        for room_num in range(1, NUM_ROOMS + 1):
            if lights[room_num] == GREEN:
                lights[room_num] = RED

    # After all guests have left, count the blue lights.
    blue_light_rooms = []
    for room_num in range(1, NUM_ROOMS + 1):
        if lights[room_num] == BLUE:
            blue_light_rooms.append(room_num)
            
    # Per the instructions, output the individual numbers that contribute to the final answer.
    # In this case, those are the room numbers of the blue lights.
    print("The following room numbers will have blue lights:")
    # Printing the list shows all the numbers involved.
    print(blue_light_rooms)
    
    # Print the final total count.
    print("\nThe final count of blue lights is:")
    print(len(blue_light_rooms))

solve_light_puzzle()
<<<32>>>