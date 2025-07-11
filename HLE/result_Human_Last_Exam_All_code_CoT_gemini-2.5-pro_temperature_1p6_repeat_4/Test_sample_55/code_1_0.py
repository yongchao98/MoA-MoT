def solve_hotel_lights_puzzle():
    """
    Solves the hotel lights puzzle by simulating the process for 100 rooms and 100 guests.
    """
    # Light states: 0 = Red, 1 = Green, 2 = Blue
    # Initial state: All 100 rooms have red lights.
    # We use a list of size 101 to use 1-based indexing for rooms 1-100.
    lights = [0] * 101

    # Loop for each guest from 1 to 100
    for guest_n in range(1, 101):
        # Guest n toggles the light n times in every nth room.
        for room_k in range(guest_n, 101, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3

        # After the guest leaves, a cat resets any green light to red.
        for room_k in range(1, 101):
            if lights[room_k] == 1:  # 1 corresponds to Green
                lights[room_k] = 0   # 0 corresponds to Red

    # After all 100 guests have passed, find which rooms have blue lights.
    blue_rooms = []
    for room_k in range(1, 101):
        if lights[room_k] == 2:  # 2 corresponds to Blue
            blue_rooms.append(room_k)

    # Output the results as described.
    print("The following rooms will have blue lights at the end:")
    # To satisfy "output each number in the final equation", we will list the room numbers.
    print(*blue_rooms, sep=", ")
    
    total_blue = len(blue_rooms)
    
    # We construct the equation string showing the sum.
    if total_blue > 0:
        equation_str = " + ".join(["1"] * total_blue)
        print(f"\nThe final count is derived from the sum of 1 for each blue room:")
        print(f"{equation_str} = {total_blue}")
    else:
        print("\nThere are no blue lights at the end.")

solve_hotel_lights_puzzle()
<<<26>>>