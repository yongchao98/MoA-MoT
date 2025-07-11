#
# Calculate the minimum storage space for a single record.
#

# --- Field Sizes in Bytes ---

# flightNumber: A standard integer and PRIMARY KEY, so it cannot be null.
# It always takes up its full space.
size_flight_number = 4  # Standard size for an integer

# ticketCost: A nullable integer. For minimum storage, we assume it's NULL.
# A NULL value occupies 0 bytes of data space.
size_ticket_cost = 0

# arrivalCity: A nullable text field. For minimum storage, we assume it's NULL.
# A NULL value occupies 0 bytes of data space.
size_arrival_city = 0

# departureCity: A nullable text field. For minimum storage, we assume it's NULL.
# A NULL value occupies 0 bytes of data space.
size_departure_city = 0

# --- Total Calculation ---
# The total minimum storage space is the sum of the minimum sizes of all fields.
# Note: This calculation does not include database-specific overhead for record headers.
total_min_storage = size_flight_number + size_ticket_cost + size_arrival_city + size_departure_city

print("Calculating the minimum storage space in bytes for a record:")
print(f"flightNumber (integer, NOT NULL): {size_flight_number} bytes")
print(f"ticketCost (integer, NULL): {size_ticket_cost} bytes")
print(f"arrivalCity (text, NULL): {size_arrival_city} bytes")
print(f"departureCity (text, NULL): {size_departure_city} bytes")
print("-" * 20)
# Print the final equation
print(f"Total = {size_flight_number} + {size_ticket_cost} + {size_arrival_city} + {size_departure_city}")
print(f"Minimum required storage space: {total_min_storage} bytes")