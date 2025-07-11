import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the FlightRecords table.
    """
    # Data type sizes in bytes
    integer_size = 4
    
    # --- Column Analysis ---
    # flightNumber is an integer and a PRIMARY KEY, so it cannot be NULL.
    flight_number_size = integer_size
    print(f"Analyzing 'flightNumber' (integer, NOT NULL): must contain data. Size: {flight_number_size} bytes.")

    # The other three fields are nullable. To minimize record size, we assume they are all NULL.
    # A NULL value does not take up space in the data portion of the record.
    ticket_cost_size = 0
    print(f"Analyzing 'ticketCost' (integer, nullable): set to NULL for minimum size. Data Size: {ticket_cost_size} bytes.")
    arrival_city_size = 0
    print(f"Analyzing 'arrivalCity' (text, nullable): set to NULL for minimum size. Data Size: {arrival_city_size} bytes.")
    departure_city_size = 0
    print(f"Analyzing 'departureCity' (text, nullable): set to NULL for minimum size. Data Size: {departure_city_size} bytes.")

    # --- Overhead Analysis ---
    # The database needs to store which columns are NULL. This is typically done with a null bitmap.
    # We have 3 nullable columns. We need 1 bit for each.
    num_nullable_columns = 3
    print(f"\nCalculating overhead for null flags...")
    print(f"Number of nullable columns: {num_nullable_columns}")

    # The space for the bitmap is the number of bits needed, rounded up to the nearest full byte.
    null_overhead_bytes = math.ceil(num_nullable_columns / 8.0)
    print(f"Overhead to store {num_nullable_columns} null flags requires {num_nullable_columns} bits, which rounds up to {null_overhead_bytes} byte(s).")

    # --- Final Calculation ---
    total_min_size = flight_number_size + null_overhead_bytes

    print("\n--- Final Calculation ---")
    print("The minimum record size is the sum of the non-nullable fields and the null flag overhead.")
    # The final print statement is formatted to show the full equation as requested.
    print(f"Final equation: {flight_number_size} (flightNumber) + {null_overhead_bytes} (null flags) = {total_min_size} bytes")

    print(f"\nMinimum storage space required: {total_min_size} bytes.")

calculate_min_record_size()
<<<5>>>