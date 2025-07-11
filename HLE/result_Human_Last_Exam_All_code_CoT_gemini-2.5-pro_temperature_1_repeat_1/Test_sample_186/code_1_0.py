import math

def solve_bomb_locations():
    """
    Calculates the highest number of bomb locations Elsa can communicate.
    """
    # Step 1: Define the initial parameters.
    video_duration_min = 5
    seconds_per_minute = 60
    # Elsa can blink or not blink each second, which is 1 bit of information per second.
    bits_per_second = 1
    map_width = 2000
    map_height = 1000

    # Step 2: Calculate the total information capacity in bits.
    total_seconds = video_duration_min * seconds_per_minute
    total_bits_capacity = total_seconds * bits_per_second

    print(f"Total time available: {video_duration_min} minutes = {total_seconds} seconds.")
    print(f"Total information Elsa can transmit: {total_seconds} seconds * {bits_per_second} bit/second = {total_bits_capacity} bits.\n")

    # Step 3: Calculate the number of bits needed to specify one location.
    total_locations = map_width * map_height
    # The number of bits needed is the ceiling of the log base 2 of the total number of locations.
    bits_per_location = math.ceil(math.log2(total_locations))

    print(f"Total possible bomb locations on the map: {map_width} * {map_height} = {total_locations}.")
    print(f"Information needed to specify one location: ceil(log2({total_locations})) = {bits_per_location} bits.\n")

    # Step 4: Calculate the maximum number of locations that can be communicated.
    # This is the total capacity divided by the cost per location.
    # We use floor because she cannot communicate a fraction of a location's coordinates.
    max_locations = math.floor(total_bits_capacity / bits_per_location)
    
    print("Final Calculation:")
    print(f"Max Locations = floor(Total Information / Bits per Location)")
    print(f"Max Locations = floor({total_bits_capacity} / {bits_per_location})")
    print(f"Max Locations = {max_locations}\n")
    print(f"The highest number of bomb locations Elsa can express is {max_locations}.")

solve_bomb_locations()
<<<14>>>