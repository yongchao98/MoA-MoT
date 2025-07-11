#
# Calculate the minimum storage space in bytes for a single record.
#

# Size of a standard integer in bytes. This field cannot be NULL.
size_flightNumber = 4

# Size for a nullable integer field set to NULL.
size_ticketCost = 0

# Size for a nullable text field set to NULL.
size_arrivalCity = 0

# Size for a nullable text field set to NULL.
size_departureCity = 0

# Calculate the total minimum size.
total_min_size = size_flightNumber + size_ticketCost + size_arrivalCity + size_departureCity

# Print the breakdown of the calculation and the final result.
print(f"Minimum storage space calculation:")
print(f"flightNumber (integer, non-nullable): {size_flightNumber} bytes")
print(f"ticketCost (integer, nullable): {size_ticketCost} bytes")
print(f"arrivalCity (text, nullable): {size_arrivalCity} bytes")
print(f"departureCity (text, nullable): {size_departureCity} bytes")
print("-" * 30)
print(f"Equation: {size_flightNumber} + {size_ticketCost} + {size_arrivalCity} + {size_departureCity} = {total_min_size}")
print(f"Total Minimum Storage: {total_min_size} bytes")
