import math

def solve_elsa_problem():
    """
    Calculates the highest number of bomb locations Elsa can express via eye blinks.
    """
    video_minutes = 5
    blink_rate_per_sec = 1
    map_width = 2000
    map_height = 1000

    # Step 1: Calculate the total time available for communication in seconds.
    total_seconds = video_minutes * 60
    print(f"Step 1: The total communication time is {video_minutes} minutes, which is {video_minutes} * 60 = {total_seconds} seconds.")

    # Step 2: Calculate the total number of bits Elsa can transmit.
    # At a rate of 1 signal per second, each second can be used to transmit one bit (blink=1, no-blink=0).
    total_available_bits = total_seconds * blink_rate_per_sec
    print(f"Step 2: With a signal rate of {blink_rate_per_sec} bit/sec, the total information Elsa can transmit is {total_seconds} * {blink_rate_per_sec} = {total_available_bits} bits.")

    # Step 3: Calculate the number of bits needed to identify one bomb location.
    # First, find the total number of unique locations on the map.
    total_unique_locations = map_width * map_height
    print(f"Step 3a: The total number of unique locations is {map_width} * {map_height} = {total_unique_locations}.")

    # To represent N unique items, we need ceil(log2(N)) bits.
    bits_per_location = math.ceil(math.log2(total_unique_locations))
    print(f"Step 3b: The bits required to specify one location is ceil(log2({total_unique_locations})) = {bits_per_location} bits.")

    # Step 4: Calculate the maximum number of locations that can be communicated.
    # This is the total bits available divided by the bits needed for each location.
    # We use integer division because she cannot communicate a fraction of a location.
    max_locations = total_available_bits // bits_per_location
    print(f"\nFinal Calculation: With {total_available_bits} total bits available and {bits_per_location} bits needed per location, the maximum number of locations is:")
    print(f"{total_available_bits} // {bits_per_location} = {max_locations}")
    
    print(f"\n<<<{max_locations}>>>")

solve_elsa_problem()