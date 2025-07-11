import math

def calculate_max_bombs():
    """
    Calculates the maximum number of bomb locations Elsa can communicate.
    """
    # --- Constants from the problem statement ---
    video_duration_minutes = 5
    blinks_per_second = 1
    map_width = 2000
    map_height = 1000

    # Step 1: Calculate the total information capacity (Total Bits)
    print("Step 1: Calculating the total bits Elsa can transmit.")
    video_duration_seconds = video_duration_minutes * 60
    total_bits_transmittable = video_duration_seconds * blinks_per_second
    print(f"Elsa can transmit for {video_duration_minutes} minutes, which is {video_duration_seconds} seconds.")
    print(f"At a rate of {blinks_per_second} bit per second, the total information capacity is: {total_bits_transmittable} bits.\n")

    # Step 2: Calculate the bits required per bomb location
    print("Step 2: Calculating the bits needed to specify one bomb location.")
    total_possible_locations = map_width * map_height
    # The number of bits needed is the ceiling of log base 2 of the number of possibilities.
    bits_per_location = math.ceil(math.log2(total_possible_locations))
    print(f"The map has {map_width}x{map_height} = {total_possible_locations} possible locations.")
    print(f"To encode one location, we need log2({total_possible_locations}) bits, which must be rounded up to {bits_per_location} bits.\n")
    
    # Step 3: Calculate the maximum number of locations
    print("Step 3: Calculating the maximum number of bomb locations.")
    # We use integer division because you can't communicate a fraction of a location.
    max_locations = total_bits_transmittable // bits_per_location
    print("This is the total transmittable bits divided by the bits needed for each location.")
    print(f"Final Equation: {total_bits_transmittable} (Total Bits) / {bits_per_location} (Bits per Location) = {max_locations} (integer division)\n")

    print(f"Therefore, the highest number of bomb locations Elsa can express is {max_locations}.")


calculate_max_bombs()
<<<14>>>