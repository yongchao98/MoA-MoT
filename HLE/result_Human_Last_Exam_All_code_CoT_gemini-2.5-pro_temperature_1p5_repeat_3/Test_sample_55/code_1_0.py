def solve_hotel_lights():
    """
    Solves the hotel lights puzzle and prints the result.
    """
    # Using 0 for Red, 1 for Green, 2 for Blue.
    # We use 101 elements to have a 1-to-1 mapping from room number to index.
    # Index 0 is ignored.
    lights = [0] * 101

    # Loop for each guest from 1 to 100
    for guest_n in range(1, 101):
        # Guest n toggles every nth room n times
        for room_number in range(guest_n, 101, guest_n):
            lights[room_number] = (lights[room_number] + guest_n) % 3

        # After each guest leaves, the cat resets any green light to red
        for room_number in range(1, 101):
            if lights[room_number] == 1:  # 1 represents Green
                lights[room_number] = 0  # 0 represents Red

    # Find and count the blue lights at the end
    blue_light_rooms = []
    for room_number in range(1, 101):
        if lights[room_number] == 2:  # 2 represents Blue
            blue_light_rooms.append(room_number)

    # Print the "equation" by showing which rooms are blue
    if not blue_light_rooms:
        print("No lights are blue at the end.")
    else:
        # Construct the output string
        output_parts = [str(room) for room in blue_light_rooms]
        print(" + ".join(output_parts) + f" = {len(blue_light_rooms)}")
        print(f"\nThe following {len(blue_light_rooms)} rooms have blue lights:")
        for room in blue_light_rooms:
            print(f"Room {room}")

    print(f"\nThe final count of blue lights is: {len(blue_light_rooms)}")


solve_hotel_lights()
<<<30>>>