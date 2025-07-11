def calculate_min_record_size():
    """
    Calculates and explains the minimum storage space for a record
    in the specified FlightRecords table.
    """

    # Size in bytes for a standard integer data type
    integer_size = 4

    # Minimum size for nullable fields is 0, as they can be NULL
    null_size = 0

    # The table has 4 columns (flightNumber, ticketCost, arrivalCity, departureCity)
    num_columns = 4
    # A null bitmap is required to track which fields are null.
    # It needs 1 bit per column. For 4 columns, we need 4 bits.
    # Storage is byte-aligned, so we need ceil(4/8) = 1 byte.
    null_bitmap_size = 1

    # flightNumber is an integer PRIMARY KEY, so it cannot be NULL.
    flight_number_size = integer_size

    # The other three fields can be NULL for the smallest record size.
    ticket_cost_size = null_size
    arrival_city_size = null_size
    departure_city_size = null_size

    # Total minimum size is the sum of the non-nullable field(s)
    # and the overhead for the null bitmap.
    total_size = (flight_number_size + ticket_cost_size +
                  arrival_city_size + departure_city_size +
                  null_bitmap_size)

    print("To find the minimum storage space, we calculate the size of the smallest possible record:")
    print(f"1. `flightNumber` (integer, non-nullable): {flight_number_size} bytes")
    print(f"2. `ticketCost` (integer, nullable): {ticket_cost_size} bytes (when NULL)")
    print(f"3. `arrivalCity` (text, nullable): {arrival_city_size} bytes (when NULL)")
    print(f"4. `departureCity` (text, nullable): {departure_city_size} bytes (when NULL)")
    print(f"5. Null Bitmap (overhead for 4 columns): {null_bitmap_size} byte")
    print("\nFinal equation for the total minimum size in bytes:")
    # Printing the final equation with each component
    print(f"{flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} + {null_bitmap_size} = {total_size}")

calculate_min_record_size()