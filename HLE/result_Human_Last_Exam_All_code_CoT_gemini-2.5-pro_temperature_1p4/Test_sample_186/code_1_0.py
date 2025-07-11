import math

def solve():
    """
    Calculates the highest number of bomb locations Elsa can express via blinking.
    """
    # Parameters from the problem
    video_duration_min = 5
    blink_rate_hz = 1  # 1 blink per second means a 1 bit/sec transmission rate
    map_width = 2000
    map_height = 1000

    # 1. Calculate total information capacity
    total_seconds = video_duration_min * 60
    # Each second can be a blink (1) or no blink (0), so the capacity is 1 bit per second.
    total_bits_capacity = total_seconds * blink_rate_hz

    # 2. Calculate bits needed per location
    # Bits for x-coordinate (0 to 1999, which is 2000 distinct values)
    bits_for_x = math.ceil(math.log2(map_width))
    # Bits for y-coordinate (0 to 999, which is 1000 distinct values)
    bits_for_y = math.ceil(math.log2(map_height))
    # Total bits to define one location
    bits_per_location = bits_for_x + bits_for_y

    # 3. Calculate the maximum number of locations
    # We can only express a whole number of locations.
    max_locations = total_bits_capacity // bits_per_location

    print(f"Total time available: {video_duration_min} minutes = {total_seconds} seconds")
    print(f"Information transmission rate: {blink_rate_hz} bit/second")
    print(f"Total information capacity: {total_seconds} seconds * {blink_rate_hz} bit/second = {total_bits_capacity} bits")
    print("-" * 20)
    print(f"Bits needed to represent the x-coordinate (up to {map_width-1}): {bits_for_x}")
    print(f"Bits needed to represent the y-coordinate (up to {map_height-1}): {bits_for_y}")
    print(f"Total bits per location: {bits_for_x} + {bits_for_y} = {bits_per_location} bits")
    print("-" * 20)
    print("Final Calculation:")
    print(f"Maximum locations = Total bits / Bits per location")
    print(f"Maximum locations = {total_bits_capacity} / {bits_per_location} = {max_locations}")

solve()
<<<14>>>