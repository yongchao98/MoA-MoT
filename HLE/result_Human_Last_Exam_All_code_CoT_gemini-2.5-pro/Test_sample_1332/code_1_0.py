# Plan:
# 1. Determine the size of the non-nullable primary key field.
# 2. Determine the number of fields that can be null.
# 3. Calculate the size of the null bitmap required to track the null status of these fields.
# 4. Sum the sizes to find the total minimum storage space for a record.

# Size in bytes for a standard integer, which cannot be null as it's a PRIMARY KEY.
flight_number_size = 4

# Number of nullable fields in the table.
# ticketCost (integer), arrivalCity (text), departureCity (text)
num_nullable_fields = 3

# To store the null status of 'num_nullable_fields', a null bitmap is used.
# We need 1 bit per nullable field. The total bytes needed is the number of bits
# divided by 8 (bits in a byte), rounded up to the nearest whole number.
# This can be calculated using integer division: (num_nullable_fields + 7) // 8
null_bitmap_size = (num_nullable_fields + 7) // 8

# The minimum total size is the sum of the non-nullable field sizes and the null bitmap size.
# The nullable fields themselves contribute 0 bytes to the data size when they are NULL.
total_min_size = flight_number_size + null_bitmap_size

print("Calculating the minimum storage space for a record:")
print(f"1. The 'flightNumber' field is a non-nullable integer and requires {flight_number_size} bytes.")
print(f"2. There are {num_nullable_fields} nullable fields. To minimize size, we assume they are all NULL.")
print(f"3. A null bitmap is needed to track these {num_nullable_fields} fields, requiring {null_bitmap_size} byte(s).")
print("\nFinal Equation:")
print(f"{flight_number_size} (for flightNumber) + {null_bitmap_size} (for null bitmap) = {total_min_size}")
