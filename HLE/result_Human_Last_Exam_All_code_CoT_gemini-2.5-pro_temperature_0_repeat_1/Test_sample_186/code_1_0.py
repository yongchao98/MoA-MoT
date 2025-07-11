import math

def solve_bomb_locations():
    """
    Calculates the maximum number of bomb locations Elsa can communicate.
    """
    # Parameters from the problem description
    video_duration_minutes = 5
    blink_rate_per_second = 1
    map_width = 2000
    map_height = 1000

    # Step 1: Calculate the total information Elsa can transmit in bits.
    # Total time in seconds = 5 minutes * 60 seconds/minute
    total_seconds = video_duration_minutes * 60
    # Total bits = total seconds * 1 bit/second (blink or no-blink)
    total_bits_available = total_seconds * blink_rate_per_second

    # Step 2: Calculate the information needed for one location in bits.
    # Bits for x-coordinate (values 0 to 1999, which is 2000 distinct values)
    # We need to find the smallest integer 'b' such that 2^b >= 2000.
    # This is calculated as ceil(log2(2000)).
    bits_for_x = math.ceil(math.log2(map_width))

    # Bits for y-coordinate (values 0 to 999, which is 1000 distinct values)
    # We need to find the smallest integer 'b' such that 2^b >= 1000.
    # This is calculated as ceil(log2(1000)).
    bits_for_y = math.ceil(math.log2(map_height))

    # Total bits per location is the sum of bits for x and y.
    bits_per_location = bits_for_x + bits_for_y

    # Step 3: Calculate the maximum number of locations.
    # This is the total available bits divided by the bits required per location.
    # We can only communicate a whole number of locations, so we take the integer floor.
    max_locations = total_bits_available // bits_per_location

    # Step 4: Print the final equation with the numbers.
    print(f"Total information Elsa can transmit: {video_duration_minutes} min * 60 s/min * {blink_rate_per_second} bit/s = {total_bits_available} bits.")
    print(f"Information needed for one location: ceil(log2({map_width})) + ceil(log2({map_height})) = {bits_for_x} + {bits_for_y} = {bits_per_location} bits.")
    print(f"Highest number of locations = floor(Total bits / Bits per location)")
    print(f"Final Equation: {total_bits_available} / {bits_per_location} = {max_locations}")

solve_bomb_locations()
<<<14>>>