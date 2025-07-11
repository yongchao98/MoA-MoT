#
# This script calculates the minimum storage space in bytes for a single record
# in the specified table structure.
#

# --- Field Size Analysis ---
# The minimum size of a record is found when all nullable fields are NULL.

# flightNumber is an integer and a PRIMARY KEY, so it cannot be NULL.
# A standard integer takes 4 bytes.
size_flightNumber = 4
field_name_flightNumber = "flightNumber (integer, PRIMARY KEY)"

# ticketCost is an integer but is nullable. Its minimum size is when it's NULL.
# A NULL value takes up 0 bytes for the data itself.
size_ticketCost = 0
field_name_ticketCost = "ticketCost (integer, NULL)"

# arrivalCity is text and is nullable. Its minimum size is when it's NULL.
size_arrivalCity = 0
field_name_arrivalCity = "arrivalCity (text, NULL)"

# departureCity is text and is nullable. Its minimum size is when it's NULL.
size_departureCity = 0
field_name_departureCity = "departureCity (text, NULL)"

# --- Calculation ---
# The total minimum size is the sum of the minimum sizes of all fields.
total_min_size = size_flightNumber + size_ticketCost + size_arrivalCity + size_departureCity

# --- Output ---
print("Calculating the minimum storage space based on data types:")
print(f"{field_name_flightNumber}: {size_flightNumber} bytes")
print(f"{field_name_ticketCost}: {size_ticketCost} bytes")
print(f"{field_name_arrivalCity}: {size_arrivalCity} bytes")
print(f"{field_name_departureCity}: {size_departureCity} bytes")
print("-" * 30)
# Final equation showing each number
print(f"Total Minimum Size = {size_flightNumber} + {size_ticketCost} + {size_arrivalCity} + {size_departureCity} = {total_min_size} bytes")
