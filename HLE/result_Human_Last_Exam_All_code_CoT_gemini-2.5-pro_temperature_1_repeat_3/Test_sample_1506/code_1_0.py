def solve():
    """
    Calculates the maximal size of the most memory-efficient data structure
    for the given prime number operations.
    """
    # The problem specifies that p and n are less than 10000.
    # This means we need to handle numbers from 0 to 9999.
    total_numbers = 10000

    # There are 8 bits in one byte.
    bits_per_byte = 8

    # The most memory-efficient data structure is a bit array where each bit
    # represents whether a number is prime or not. This allows for O(1) lookups.
    #
    # To calculate the size of this bit array in bytes, we need to find out
    # how many bytes are required to store `total_numbers` bits.
    # We use ceiling division, as we can't have a fraction of a byte.
    # For integers, ceil(a / b) is equivalent to (a + b - 1) // b.
    # Since 10000 is perfectly divisible by 8, simple division also works.
    size_in_bytes = total_numbers // bits_per_byte

    print("Step 1: Choose the most memory-efficient data structure.")
    print("A bit array is the best choice for O(1) primality tests and memory efficiency.")
    print("\nStep 2: Calculate the size of this data structure.")
    print(f"We need to store primality information for {total_numbers} numbers (0 to 9999).")
    print("Each number requires 1 bit.")
    print(f"Total bits needed = {total_numbers} bits.")
    print(f"There are {bits_per_byte} bits in a byte.")
    print("\nThe final equation for the size in bytes is:")
    print(f"{total_numbers} (bits) / {bits_per_byte} (bits per byte) = {size_in_bytes} (bytes)")

solve()