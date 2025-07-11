import math

def solve_bomb_locations():
    """
    Calculates the maximum number of bomb locations Elsa can communicate.
    """
    # Step 1: Calculate the total information capacity (in bits)
    video_duration_min = 5
    blink_rate_per_sec = 1
    
    total_seconds = video_duration_min * 60
    # Each second can represent one bit of information (blink or no blink)
    total_bits_available = total_seconds * blink_rate_per_sec
    
    # Step 2: Calculate the information needed per location (in bits)
    map_width = 2000
    map_height = 1000
    
    # Bits needed to represent the x-coordinate (0 to 1999)
    bits_for_x = math.ceil(math.log2(map_width))
    
    # Bits needed to represent the y-coordinate (0 to 999)
    bits_for_y = math.ceil(math.log2(map_height))
    
    # Total bits for one location (x, y)
    bits_per_location = bits_for_x + bits_for_y
    
    # Step 3: Calculate the maximum number of locations
    # We use integer division as Elsa can only communicate whole locations.
    max_locations = total_bits_available // bits_per_location
    
    # Print the step-by-step calculation
    print(f"The video is {video_duration_min} minutes long, which is {total_seconds} seconds.")
    print(f"With a maximum of {blink_rate_per_sec} blink per second, Elsa can transmit a total of {total_bits_available} bits of information.")
    print("-" * 20)
    print(f"The map is {map_width}x{map_height}.")
    print(f"To encode an x-coordinate (up to {map_width-1}), we need ceil(log2({map_width})) = {bits_for_x} bits.")
    print(f"To encode a y-coordinate (up to {map_height-1}), we need ceil(log2({map_height})) = {bits_for_y} bits.")
    print(f"Total bits to specify one bomb location = {bits_for_x} + {bits_for_y} = {bits_per_location} bits.")
    print("-" * 20)
    print("To find the highest number of locations, we divide the total bits available by the bits per location:")
    print(f"Equation: {total_bits_available} // {bits_per_location} = {max_locations}")
    print(f"\nTherefore, the highest number of bomb locations Elsa can express is {max_locations}.")

solve_bomb_locations()
<<<14>>>