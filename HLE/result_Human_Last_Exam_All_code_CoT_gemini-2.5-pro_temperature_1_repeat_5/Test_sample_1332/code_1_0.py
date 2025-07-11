# Step 1: Define the size in bytes for the data types involved.
# A standard integer typically requires 4 bytes.
integer_size = 4

# Step 2: Account for the data storage of the non-nullable field.
# The 'flightNumber' is an integer and a PRIMARY KEY, so it cannot be NULL.
# It will always consume its full data size.
flight_number_size = integer_size
print(f"The 'flightNumber' field is a non-nullable integer and requires {flight_number_size} bytes.")

# Step 3: Account for the data storage of the nullable fields for a minimum record.
# To get the minimum record size, we assume all nullable fields are NULL.
# NULL values do not consume data space for their respective columns.
ticket_cost_min_size = 0
arrival_city_min_size = 0
departure_city_min_size = 0
print("For minimum record size, 'ticketCost', 'arrivalCity', and 'departureCity' are NULL, contributing 0 data bytes each.")

# Step 4: Calculate the overhead required to track which columns are NULL.
# A null bitmap is used, with one bit per column.
# The table has 4 columns, so we need 4 bits.
# Storage is allocated in whole bytes, so 1 byte (8 bits) is sufficient.
null_bitmap_size = 1
print(f"A null bitmap is required to track NULLs, adding {null_bitmap_size} byte of overhead.")

# Step 5: Calculate the total minimum storage space.
# This is the sum of the non-nullable data and the null bitmap overhead.
total_min_size = flight_number_size + null_bitmap_size

# Final Output
print("\nTo find the minimum storage space, we sum the size of the mandatory data and the null tracking overhead.")
print(f"Final Calculation: {flight_number_size} bytes (flightNumber) + {null_bitmap_size} byte (null bitmap) = {total_min_size} bytes")
print(f"The minimum storage space required for a single record is {total_min_size} bytes.")
<<<5>>>