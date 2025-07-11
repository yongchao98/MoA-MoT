def solve_hotel_lights_puzzle():
    """
    Simulates the hotel lights problem to find the number of blue lights.
    """
    # Step 1: Initialization
    # 0: red, 1: green, 2: blue
    # Using a list of size 101 for 1-based indexing (rooms 1-100)
    lights = [0] * 101

    # Step 2 & 3: Guest and Cat Simulation
    # Loop through guests from 1 to 100
    for guest_n in range(1, 101):
        # Guest n toggles every nth room n times
        for room_number in range(guest_n, 101, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3

        # The cat resets any green light to red
        for room_number in range(1, 101):
            if lights[room_number] == 1:  # if green
                lights[room_number] = 0  # reset to red

    # Step 4 & 5: Final Count and Output
    # Create a list of indicators (1 for blue, 0 for not blue) for the final equation.
    blue_light_indicators = []
    for i in range(1, 101):
        if lights[i] == 2:
            blue_light_indicators.append("1")
        else:
            blue_light_indicators.append("0")

    # Calculate the total number of blue lights
    blue_count = blue_light_indicators.count("1")

    # Print the "equation" as requested, showing a 1 for each blue light
    print("This equation represents each room from 1 to 100. A '1' means the light is blue, a '0' means it is not.")
    equation_str = " + ".join(blue_light_indicators)
    print(f"{equation_str} = {blue_count}")
    print(f"\nThere are a total of {blue_count} blue lights at the end.")

solve_hotel_lights_puzzle()
<<<33>>>