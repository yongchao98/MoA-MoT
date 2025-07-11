# Plan:
# 1. Define the size in bytes for the non-nullable 'flightNumber' field.
# 2. Define the size in bytes for the null bitmap needed to track the 3 nullable fields.
# 3. Sum these values to get the total minimum storage space.
# 4. Print the calculation and the final result.

# Size of a standard integer in bytes.
# The 'flightNumber' is a PRIMARY KEY, so it cannot be NULL and will always occupy this space.
flight_number_size = 4

# Number of nullable fields in the table.
# These are: ticketCost, arrivalCity, departureCity.
num_nullable_fields = 3

# To store the NULL/NOT NULL status for 3 fields, we need 3 bits.
# Since storage is byte-aligned, these bits are stored in a single byte.
null_bitmap_size = 1

# The minimum storage for the record is the sum of the non-nullable field's size
# and the size of the null bitmap. The nullable fields themselves contribute 0 bytes
# to the data size when they are NULL.
total_minimum_size = flight_number_size + null_bitmap_size

print("Calculating the minimum storage space for a record:")
print(f"Size of non-nullable 'flightNumber' (integer): {flight_number_size} bytes")
print(f"Size of null bitmap (to track {num_nullable_fields} nullable fields): {null_bitmap_size} byte")
print(f"Size of nullable fields ('ticketCost', 'arrivalCity', 'departureCity') when NULL: 0 bytes")
print("\nFinal Equation:")
print(f"{flight_number_size} (flightNumber) + {null_bitmap_size} (null bitmap) = {total_minimum_size}")

print(f"\nThe minimum storage space required is {total_minimum_size} bytes.")