import math

# --- Plan ---
# 1. Define the size of the non-nullable field (flightNumber).
# 2. Define the total number of columns to calculate the NULL bitmap size.
# 3. Assume all nullable fields are NULL to get the minimum size.
# 4. The minimum size will be the size of the non-nullable data plus the size of the NULL bitmap.

# Size of a standard integer in bytes
flight_number_size = 4

# Number of columns in the table
total_columns = 4

# A NULL bitmap is used to track which columns are NULL.
# It needs 1 bit per column.
bits_for_null_bitmap = total_columns

# Storage is byte-addressable. Calculate bytes needed to store the bits.
# We use ceiling division: math.ceil(bits / 8) or integer division (bits + 7) // 8
null_bitmap_size_bytes = (bits_for_null_bitmap + 7) // 8

# For the smallest possible record, all nullable fields are NULL.
# Their data size is 0. The only data stored is for the non-nullable PRIMARY KEY.
# Total minimum size = (size of fixed non-nullable data) + (size of null bitmap)
total_min_size = flight_number_size + null_bitmap_size_bytes

# --- Output ---
print("To find the minimum storage space for a record, we calculate the size of mandatory fields plus any structural overhead.")
print("1. The 'flightNumber' is a non-nullable integer, so it always requires 4 bytes.")
print(f"Size of 'flightNumber' field: {flight_number_size} bytes")

print("\n2. The other three fields can be NULL. To store this NULL information, a NULL bitmap is used.")
print(f"The table has {total_columns} columns, requiring a {bits_for_null_bitmap}-bit bitmap. This occupies {null_bitmap_size_bytes} byte(s).")
print(f"Overhead for NULL bitmap: {null_bitmap_size_bytes} byte(s)")

print("\n3. The total minimum storage space is the sum of the mandatory data and the NULL bitmap.")
print(f"Final Equation: {flight_number_size} (for flightNumber) + {null_bitmap_size_bytes} (for NULL bitmap) = {total_min_size} bytes")
