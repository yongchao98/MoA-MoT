import math

def solve_hotel_lights():
    """
    Solves the hotel lights puzzle and prints the number of blue lights.
    
    States are represented numerically:
    0: Red
    1: Green
    2: Blue
    """

    # We will determine the final state of each room's light one by one.
    blue_lights_count = 0

    # Iterate through each room from 1 to 100.
    for room_number in range(1, 101):
        # The light in room 'room_number' is only affected by guests whose number
        # is a divisor of the room number.
        
        # Step 1: Find all divisors of the current room_number.
        divisors = []
        for i in range(1, int(math.sqrt(room_number)) + 1):
            if room_number % i == 0:
                divisors.append(i)
                if i * i != room_number:
                    divisors.append(room_number // i)
        divisors.sort() # Guests arrive in order, so we process divisors in order.

        # Step 2: Simulate the light changes for this single room.
        # Initial state is Red (0).
        current_state = 0 
        
        for guest_number in divisors:
            # Guest 'guest_number' toggles the light 'guest_number' times.
            # A toggle moves to the next color in the cycle R->G->B->R.
            # Toggling 'n' times is equivalent to adding 'n' to the state (mod 3).
            current_state = (current_state + guest_number) % 3
            
            # After each guest leaves, a cat resets any Green (1) light to Red (0).
            if current_state == 1:
                current_state = 0
        
        # Step 3: Check the final state of the room's light.
        if current_state == 2: # If the light is Blue
            blue_lights_count += 1
    
    print(f"The number of lights that will be blue at the end is: {blue_lights_count}")

solve_hotel_lights()
<<<30>>>