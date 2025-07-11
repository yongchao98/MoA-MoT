import math

def calculate_min_storage():
    """
    Calculates and explains the minimum storage space for a record based on the
    provided SQL table structure.

    The calculation makes the following standard assumptions:
    - An 'integer' data type occupies 4 bytes.
    - A 'PRIMARY KEY' field cannot be NULL.
    - A NULL value occupies 0 bytes of data storage.
    - A null bitmap is used to track NULL values, requiring 1 bit per
      nullable column, rounded up to the nearest byte.
    """
    # Size in bytes for a standard integer.
    integer_size = 4

    # The 'flightNumber' is an integer and cannot be NULL.
    flight_number_size = integer_size

    # For minimum record size, nullable fields are assumed to be NULL.
    # Their data storage size is 0 bytes.
    ticket_cost_size = 0
    arrival_city_size = 0
    departure_city_size = 0

    # There are 3 nullable columns (ticketCost, arrivalCity, departureCity).
    # The overhead to track these NULLs (e.g., via a null bitmap) is
    # calculated by taking 1 bit per column and rounding up to the nearest byte.
    num_nullable_columns = 3
    null_overhead_bytes = math.ceil(num_nullable_columns / 8.0)

    # The total minimum size is the sum of the non-nullable field sizes
    # and the null-tracking overhead.
    total_min_size = (flight_number_size +
                      ticket_cost_size +
                      arrival_city_size +
                      departure_city_size +
                      null_overhead_bytes)

    print("Calculation of Minimum Record Storage:")
    print("---------------------------------------")
    print(f"Size of 'flightNumber' (non-nullable INTEGER): {flight_number_size} bytes")
    print(f"Size of 'ticketCost' (when NULL): {ticket_cost_size} bytes")
    print(f"Size of 'arrivalCity' (when NULL): {arrival_city_size} bytes")
    print(f"Size of 'departureCity' (when NULL): {departure_city_size} bytes")
    print(f"Overhead for tracking 3 NULL fields (null bitmap): {int(null_overhead_bytes)} byte")
    print("\nFinal Equation (sum of all parts):")
    print(f"{flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} + {int(null_overhead_bytes)} = {int(total_min_size)}")

calculate_min_storage()
<<<5>>>