import math

# This script calculates the size of the most memory-efficient data structure
# for the described prime number operations.

# The design uses a bit array to store the primality of numbers up to the limit.
# The maximum number, p, is less than 10000, so we need to handle numbers 0 through 9999.
limit = 10000

# Each number's primality (prime or not) can be stored in a single bit.
bits_per_number = 1
total_bits = limit * bits_per_number

# There are 8 bits in a standard byte.
bits_in_a_byte = 8

# To find the size in bytes, we divide the total number of bits required by the
# number of bits in a byte. We use math.ceil to ensure we have enough bytes
# to store all the bits, as memory is allocated in whole bytes.
size_in_bytes = math.ceil(total_bits / bits_in_a_byte)

print("--- Data Structure Size Calculation ---")
print(f"The design stores primality for N = {limit} numbers (from 0 to 9999).")
print(f"To be memory-efficient, each number's primality is stored using {bits_per_number} bit.")
print(f"Total bits required = {total_bits}")
print(f"Number of bits in one byte = {bits_in_a_byte}")
print("")
print("The final equation to calculate the size in bytes is:")
print(f"   ceil({total_bits} / {bits_in_a_byte})")
print("")
print(f"The maximal size of this data structure is {int(size_in_bytes)} bytes.")
