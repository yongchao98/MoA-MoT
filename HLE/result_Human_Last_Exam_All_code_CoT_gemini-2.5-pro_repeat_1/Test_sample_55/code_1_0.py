def solve_hotel_lights_puzzle():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    NUM_ROOMS = 100
    
    # Light states: 0 for Red, 1 for Green, 2 for Blue
    # Initialize all 100 rooms to Red. Use a list of size 101 to map room number to index directly.
    lights = [0] * (NUM_ROOMS + 1)

    # Loop for each guest from 1 to 100
    for guest_n in range(1, NUM_ROOMS + 1):
        # Guest n toggles the light in every nth room n times
        for room_number in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3
            
        # After each guest leaves, the cat resets any green light to red
        for i in range(1, NUM_ROOMS + 1):
            if lights[i] == 1:  # If the light is Green
                lights[i] = 0   # Reset to Red

    # After all guests have left, find the rooms with blue lights
    blue_light_rooms = []
    for i in range(1, NUM_ROOMS + 1):
        if lights[i] == 2:  # If the light is Blue
            blue_light_rooms.append(i)

    # As requested, output each number in the final equation (the room numbers)
    print("The rooms with blue lights are:")
    # Printing the numbers themselves, not as a list representation
    for room_num in blue_light_rooms:
        print(room_num, end=' ')
    print("\n") # Newline for better formatting

    # Output the final count
    final_count = len(blue_light_rooms)
    print(f"The total number of blue lights at the end is: {final_count}")

solve_hotel_lights_puzzle()
<<<32>>>