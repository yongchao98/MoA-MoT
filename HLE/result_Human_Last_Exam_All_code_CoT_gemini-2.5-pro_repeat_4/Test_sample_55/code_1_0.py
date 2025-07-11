def solve_hotel_lights_puzzle():
    """
    Solves the hotel lights puzzle by simulating the actions of 100 guests and a cat.
    """
    # Use integers to represent the state of the lights for easier cycling.
    # 0: Red, 1: Green, 2: Blue
    RED, GREEN, BLUE = 0, 1, 2
    
    num_rooms = 100
    
    # Initialize all 100 lights to Red.
    # We use a list of size 101 to allow for 1-based indexing (rooms 1 to 100).
    lights = [RED] * (num_rooms + 1)
    
    # Loop through each guest from 1 to 100.
    for guest_n in range(1, num_rooms + 1):
        
        # Guest 'n' toggles the light in every 'n'th room.
        for room_number in range(guest_n, num_rooms + 1, guest_n):
            # Toggling 'n' times is equivalent to adding 'n' to the current state
            # and taking the result modulo 3 to handle the color cycle.
            lights[room_number] = (lights[room_number] + guest_n) % 3
            
        # After each guest leaves, the cat resets any Green light to Red.
        for i in range(1, num_rooms + 1):
            if lights[i] == GREEN:
                lights[i] = RED
                
    # After all guests have passed, find which rooms have a Blue light.
    blue_light_rooms = []
    for i in range(1, num_rooms + 1):
        if lights[i] == BLUE:
            blue_light_rooms.append(i)
            
    # As requested, output the numbers that make up the final answer.
    # First, we print the room numbers themselves.
    print("The room numbers of the lights that are blue are:")
    # The str() conversion is to ensure a clean print of the list of numbers.
    print(str(blue_light_rooms))
    
    # Then, we show the "equation" for the final count.
    print("\nThe equation for the final count is a sum of '1' for each blue room:")
    if not blue_light_rooms:
        # Handle the case where no lights are blue.
        print(f"0 = 0")
    else:
        equation_str = " + ".join(["1"] * len(blue_light_rooms))
        print(f"{equation_str} = {len(blue_light_rooms)}")

solve_hotel_lights_puzzle()
<<<20>>>