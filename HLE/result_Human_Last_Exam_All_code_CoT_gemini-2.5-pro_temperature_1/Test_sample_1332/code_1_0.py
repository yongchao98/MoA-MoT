# Step 1: Define the storage size for the non-nullable primary key.
# A standard integer requires 4 bytes. As a PRIMARY KEY, it cannot be NULL.
flightNumber_size = 4

# Step 2: Define the minimum storage size for the nullable fields.
# For a record to be as small as possible, all nullable fields will be set to NULL.
# A NULL value itself does not occupy data space for the column.
ticketCost_min_size = 0
arrivalCity_min_size = 0
departureCity_min_size = 0

# Step 3: Calculate the total minimum storage space by summing the sizes.
total_min_size = flightNumber_size + ticketCost_min_size + arrivalCity_min_size + departureCity_min_size

# Step 4: Print the explanation and the final calculation.
print("To find the minimum storage space, we sum the minimum size of each field:")
print(f"flightNumber (integer, PRIMARY KEY): {flightNumber_size} bytes (cannot be NULL)")
print(f"ticketCost (integer, nullable): {ticketCost_min_size} bytes (when NULL)")
print(f"arrivalCity (text, nullable): {arrivalCity_min_size} bytes (when NULL)")
print(f"departureCity (text, nullable): {departureCity_min_size} bytes (when NULL)")
print("-" * 20)
print("Total Minimum Storage Calculation:")
print(f"{flightNumber_size} + {ticketCost_min_size} + {arrivalCity_min_size} + {departureCity_min_size} = {total_min_size}")
print(f"\nThe minimum storage space required for a single record is {total_min_size} bytes.")
