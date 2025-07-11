# Step 1: Define the size in bytes for the data types in the smallest possible record.
# A standard integer takes 4 bytes. As a PRIMARY KEY, 'flightNumber' cannot be NULL.
flightNumber_size = 4

# The other fields can be NULL. For the minimum record size, we assume they are NULL.
# The storage space for a NULL value in the data portion of a record is 0 bytes.
null_value_size = 0
ticketCost_size = null_value_size
arrivalCity_size = null_value_size
departureCity_size = null_value_size

# Step 2: Calculate the total minimum storage space for one record.
total_min_size = flightNumber_size + ticketCost_size + arrivalCity_size + departureCity_size

# Step 3: Print the equation and the final result.
print(f"The minimum storage space is calculated by summing the sizes of each field in its smallest possible state.")
print(f"Minimum size of flightNumber (integer, NOT NULL): {flightNumber_size} bytes")
print(f"Minimum size of ticketCost (integer, nullable): {ticketCost_size} bytes")
print(f"Minimum size of arrivalCity (text, nullable): {arrivalCity_size} bytes")
print(f"Minimum size of departureCity (text, nullable): {departureCity_size} bytes")
print(f"\nFinal Equation:")
print(f"{flightNumber_size} + {ticketCost_size} + {arrivalCity_size} + {departureCity_size} = {total_min_size}")