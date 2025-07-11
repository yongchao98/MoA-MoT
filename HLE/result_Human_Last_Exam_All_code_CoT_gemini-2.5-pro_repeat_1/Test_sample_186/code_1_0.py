import math

def solve_bomb_locations():
    """
    Calculates the highest number of bomb locations Elsa can express.
    """
    # Parameters from the problem description
    map_width = 2000
    map_height = 1000
    video_duration_min = 5
    blinks_per_second = 1  # This represents the information rate in bits/second

    # Step 1: Calculate the total number of bits Elsa can transmit.
    # Total time in seconds = 5 minutes * 60 seconds/minute
    total_seconds = video_duration_min * 60
    # Total bits available = total time * information rate
    total_bits_available = total_seconds * blinks_per_second

    # Step 2: Calculate the number of bits needed to encode one bomb location.
    # A location is (x,y) where 0 <= x < 2000 and 0 <= y < 1000.
    # We need to find the number of bits to represent 2000 distinct values for x
    # and 1000 distinct values for y.
    # Bits needed for a number 'N' is ceil(log2(N)).
    bits_for_x = math.ceil(math.log2(map_width))
    bits_for_y = math.ceil(math.log2(map_height))
    
    # Total bits per location is the sum of bits for x and y.
    bits_per_location = bits_for_x + bits_for_y

    # Step 3: Calculate the maximum number of locations.
    # This is the integer division of total available bits by bits per location.
    max_locations = total_bits_available // bits_per_location

    # Print the steps of the calculation as requested.
    print(f"1. Total communication time: {video_duration_min} min * 60 s/min = {total_seconds} seconds.")
    print(f"2. Total information capacity: {total_seconds} seconds * {blinks_per_second} bit/second = {total_bits_available} bits.")
    print(f"3. Information per location:")
    print(f"   - Bits for x-coordinate (range 0-{map_width-1}): ceil(log2({map_width})) = {bits_for_x} bits.")
    print(f"   - Bits for y-coordinate (range 0-{map_height-1}): ceil(log2({map_height})) = {bits_for_y} bits.")
    print(f"   - Total bits per location = {bits_for_x} + {bits_for_y} = {bits_per_location} bits.")
    print(f"\n4. Final Calculation:")
    print(f"   Highest number of locations = floor(Total bits available / Bits per location)")
    print(f"   Highest number of locations = floor({total_bits_available} / {bits_per_location}) = {max_locations}")

solve_bomb_locations()
<<<14>>>