import math

# 1. Define the size of the mandatory, non-nullable field.
# 'flightNumber' is an integer PRIMARY KEY and cannot be null. A standard integer uses 4 bytes.
size_of_flightNumber = 4

# 2. Determine the number of nullable fields.
# 'ticketCost', 'arrivalCity', and 'departureCity' can be NULL.
# For the minimum record size, their data storage contribution is 0 bytes.
num_nullable_fields = 3

# 3. Calculate the storage required for the null-tracking mechanism (null bitmap).
# We need one bit for each nullable field.
bits_for_null_bitmap = num_nullable_fields
# Storage is allocated in whole bytes, so we find the ceiling of bits / 8.
size_of_null_bitmap = math.ceil(bits_for_null_bitmap / 8)

# 4. Calculate the total minimum storage size.
# This is the sum of the non-nullable data and the null bitmap.
total_minimum_size = size_of_flightNumber + size_of_null_bitmap

# 5. Print the final result and the equation.
print("This script calculates the minimum storage space for a record by summing the size of non-nullable fields and the overhead required to track null values.")
print("\n--- Calculation Breakdown ---")
print(f"Size of mandatory 'flightNumber' (integer): {size_of_flightNumber} bytes")
print(f"Number of nullable fields: {num_nullable_fields}")
print(f"Storage for null bitmap ({num_nullable_fields} bits required): {size_of_null_bitmap} byte(s)")
print("\n--- Final Equation ---")
print(f"Total Minimum Storage = (Size of 'flightNumber') + (Size of Null Bitmap)")
# Final print of the equation with the numbers
print(f"{int(total_minimum_size)} = {size_of_flightNumber} + {int(size_of_null_bitmap)}")
