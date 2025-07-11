# Step 1: Define the size for the non-nullable primary key.
# A standard integer requires 4 bytes. As a PRIMARY KEY, flightNumber cannot be null.
flightNumber_size = 4

# Step 2: Define the minimum size for the nullable fields.
# To achieve the smallest record size, all nullable fields are assumed to be NULL.
# A NULL value does not occupy space in the data portion of a row.
ticketCost_min_size = 0
arrivalCity_min_size = 0
departureCity_min_size = 0

# Step 3: Calculate the total minimum storage space by summing the sizes.
total_min_storage = flightNumber_size + ticketCost_min_size + arrivalCity_min_size + departureCity_min_size

# Step 4: Print the final equation and the result.
print(f"The minimum storage is calculated by summing the fixed size of the non-nullable primary key and assuming all other nullable fields are NULL.")
print(f"Minimum Storage = (size of flightNumber) + (size of ticketCost) + (size of arrivalCity) + (size of departureCity)")
print(f"Calculation: {flightNumber_size} + {ticketCost_min_size} + {arrivalCity_min_size} + {departureCity_min_size} = {total_min_storage} bytes")