def solve_hotel_lights_riddle():
    """
    This function simulates the hotel light riddle to find the number of blue lights.
    - There are 100 rooms, lights can be Red, Green, or Blue.
    - We map states: 0=Red, 1=Green, 2=Blue.
    - Initially, all 100 lights are Red (0).
    - For each guest `n` from 1 to 100:
        - Guest `n` toggles the light in every nth room `n` times.
          This means state becomes (state + n) % 3.
        - A cat resets any Green (1) light to Red (0).
    - The function counts how many lights are Blue (2) at the end.
    """
    num_rooms = 100
    
    # Initialize lights for rooms 1 to 100. Index 0 is unused.
    # All lights start as Red (0).
    lights = [0] * (num_rooms + 1)

    # Loop through each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        
        # Guest `n` visits every nth room and toggles the light `n` times
        for room_number in range(guest_n, num_rooms + 1, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3

        # After the guest leaves, the cat resets any green light to red
        for room_number in range(1, num_rooms + 1):
            if lights[room_number] == 1:  # If light is Green
                lights[room_number] = 0   # Reset to Red

    # Count the number of blue lights at the end
    blue_count = 0
    for room_number in range(1, num_rooms + 1):
        if lights[room_number] == 2:  # If light is Blue
            blue_count += 1
    
    # The final equation is the count of blue lights.
    print(f"Final count of blue lights = {blue_count}")


solve_hotel_lights_riddle()
<<<30>>>