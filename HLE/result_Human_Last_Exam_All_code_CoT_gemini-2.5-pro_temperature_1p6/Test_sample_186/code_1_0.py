import math

def solve_bomb_locations():
    """
    Calculates the highest number of bomb locations Elsa can express.
    """
    # 1. Calculate the total information capacity (in bits).
    # The video is 5 minutes long.
    video_duration_min = 5
    # She can transmit information at a rate of 1 blink per second, which is equivalent to 1 bit per second.
    info_rate_bps = 1

    total_seconds = video_duration_min * 60
    total_bits_available = total_seconds * info_rate_bps

    print(f"Total communication time: {video_duration_min} minutes = {total_seconds} seconds.")
    print(f"Total information capacity: {total_seconds} seconds * {info_rate_bps} bit/second = {int(total_bits_available)} bits.\n")

    # 2. Calculate the information needed per location (in bits).
    # The map resolution is 2000x1000.
    map_x = 2000
    map_y = 1000

    # Bits needed for X coordinate (from 0 to 1999)
    bits_for_x = math.ceil(math.log2(map_x))
    # Bits needed for Y coordinate (from 0 to 999)
    bits_for_y = math.ceil(math.log2(map_y))
    
    total_bits_per_location = bits_for_x + bits_for_y

    print(f"To specify an X coordinate up to {map_x}, we need ceil(log2({map_x})) = {bits_for_x} bits.")
    print(f"To specify a Y coordinate up to {map_y}, we need ceil(log2({map_y})) = {bits_for_y} bits.")
    print(f"Total bits per location = {bits_for_x} + {bits_for_y} = {total_bits_per_location} bits.\n")

    # 3. Determine the maximum number of locations.
    # This is the total bits available divided by the bits needed for one location.
    # We use floor because she can only communicate a whole number of locations.
    max_locations = math.floor(total_bits_available / total_bits_per_location)
    
    print("The maximum number of locations is the total bits available divided by the bits per location.")
    # The final print statement matches the required format.
    print(f"Maximum locations = floor({int(total_bits_available)} / {total_bits_per_location}) = {max_locations}")
    
    return max_locations

# Execute the function to get the final answer.
final_answer = solve_bomb_locations()
# The final answer is wrapped according to the format requirement.
print(f"\n<<<{final_answer}>>>")