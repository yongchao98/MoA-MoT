def solve_hotel_lights_puzzle():
    """
    Solves the hotel lights puzzle by simulating the actions of guests and a cat.
    """
    # Define states for clarity: 0=Red, 1=Green, 2=Blue
    RED, GREEN, BLUE = 0, 1, 2
    
    # Initialize 100 rooms (plus an unused 0th room for 1-based indexing)
    # All lights start as RED.
    lights = [RED] * 101

    # Loop through each guest from 1 to 100
    for guest_n in range(1, 101):
        # Guest n toggles every nth room n times
        for room_number in range(guest_n, 101, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3
            
        # After each guest, the cat resets any green lights to red
        for room_number in range(1, 101):
            if lights[room_number] == GREEN:
                lights[room_number] = RED

    # Find all rooms that are blue at the end
    blue_rooms = []
    for room_number in range(1, 101):
        if lights[room_number] == BLUE:
            blue_rooms.append(str(room_number))
            
    # Print the final result as an equation
    if not blue_rooms:
        print("No lights ended up being blue.")
        print("0 = 0")
    else:
        # Create a string of "1 + 1 + ..." for each blue light
        sum_of_ones = " + ".join(["1"] * len(blue_rooms))
        total_blue = len(blue_rooms)
        
        print(f"The puzzle is solved. The final state of the lights has been calculated.")
        print(f"The rooms with blue lights are: {', '.join(blue_rooms)}.")
        print(f"Counting each of these {total_blue} rooms gives the following equation:")
        print(f"{sum_of_ones} = {total_blue}")

# Run the simulation and print the result
solve_hotel_lights_puzzle()
<<<23>>>