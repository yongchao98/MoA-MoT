#
# Calculate the minimum storage space in bytes for a single record.
#

# The 'flightNumber' is an integer PRIMARY KEY. It cannot be NULL.
# A standard integer requires 4 bytes of storage.
flight_number_size = 4

# The 'ticketCost' field is an integer but can be NULL.
# For the smallest possible record, we assume it is NULL.
# A NULL value takes up 0 bytes of data space.
ticket_cost_size = 0

# The 'arrivalCity' field is text and can be NULL.
# For minimum storage, we assume it is NULL.
# A NULL value takes up 0 bytes of data space.
arrival_city_size = 0

# The 'departureCity' field is text and can also be NULL.
# For minimum storage, we assume it is NULL.
# A NULL value takes up 0 bytes of data space.
departure_city_size = 0

# The total minimum size is the sum of the sizes of all fields.
total_size = flight_number_size + ticket_cost_size + arrival_city_size + departure_city_size

# Print the breakdown of the calculation and the final result.
print("Calculating the minimum storage space:")
print(f"flightNumber (integer, PRIMARY KEY): {flight_number_size} bytes")
print(f"ticketCost (nullable integer, as NULL): {ticket_cost_size} bytes")
print(f"arrivalCity (nullable text, as NULL): {arrival_city_size} bytes")
print(f"departureCity (nullable text, as NULL): {departure_city_size} bytes")
print("-" * 20)
print(f"Total Minimum Size = {flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} = {total_size} bytes")
