import math

# Plan:
# 1. Define the size of each data type. An 'integer' typically takes 4 bytes.
# 2. For the smallest record, all nullable fields ('ticketCost', 'arrivalCity', 'departureCity') are NULL.
#    NULL values do not consume space in the data area of the record.
# 3. The 'flightNumber' is a PRIMARY KEY, so it cannot be NULL. It will always consume its data type's space.
# 4. Calculate the size of the required data. This will be the size of 'flightNumber'.
# 5. Calculate the overhead for tracking NULLs. A database uses a "null bitmap" for this.
#    There are 4 columns in total, so we need 4 bits. Since storage is byte-aligned,
#    we calculate math.ceil(4 / 8) to find the bytes needed for the bitmap.
# 6. Sum the data size and the null bitmap size for the final answer.

# --- Calculation ---

# Size of the mandatory integer field
size_flightNumber = 4  # bytes

# Number of columns in the table
num_columns = 4

# Size of the null bitmap
# We need 1 bit for each of the 4 columns to indicate if it's null or not.
# Storage is allocated in whole bytes.
size_null_bitmap = math.ceil(num_columns / 8) # This equals 1 byte

# The total minimum size is the sum of the data that must be stored
# and the overhead to track which columns are null.
total_min_size = size_flightNumber + size_null_bitmap

print("Calculating the minimum storage space for a record:")
print(f"1. Non-nullable data 'flightNumber' (integer): {size_flightNumber} bytes")
print(f"2. Null bitmap for all {num_columns} columns (math.ceil({num_columns}/8)): {size_null_bitmap} byte")
print("\nFinal equation:")
print(f"{size_flightNumber} + {size_null_bitmap} = {total_min_size}")