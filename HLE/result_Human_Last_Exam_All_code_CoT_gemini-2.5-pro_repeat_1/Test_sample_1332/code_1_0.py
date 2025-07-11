def calculate_minimum_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the specified table.
    """
    # Step 1: Define the size of data types in bytes.
    # A standard SQL INTEGER typically requires 4 bytes.
    integer_size = 4

    # Step 2: Analyze each column for its minimum possible size.

    # 'flightNumber' is a non-nullable integer (PRIMARY KEY).
    # Its size is fixed.
    flight_number_min_size = integer_size
    print(f"Size of 'flightNumber' (integer, non-nullable): {flight_number_min_size} bytes")

    # For nullable columns, the smallest representation is NULL.
    # A NULL value itself takes up 0 bytes in the data section of the record.
    ticket_cost_min_size = 0
    print(f"Minimum size of 'ticketCost' (when NULL): {ticket_cost_min_size} bytes")

    arrival_city_min_size = 0
    print(f"Minimum size of 'arrivalCity' (when NULL): {arrival_city_min_size} bytes")

    departure_city_min_size = 0
    print(f"Minimum size of 'departureCity' (when NULL): {departure_city_min_size} bytes")

    # Step 3: Calculate the overhead for storing NULL information.
    # Databases use a null bitmap to track NULLs. It needs 1 bit per nullable column.
    # We have 3 nullable columns (ticketCost, arrivalCity, departureCity).
    # This requires 3 bits. Since storage is byte-aligned, this takes 1 byte.
    num_nullable_columns = 3
    null_bitmap_size = 1
    print(f"Overhead size for null bitmap ({num_nullable_columns} nullable columns): {null_bitmap_size} byte")

    # Step 4: Sum the sizes to get the total minimum record size.
    total_min_size = (flight_number_min_size +
                      ticket_cost_min_size +
                      arrival_city_min_size +
                      departure_city_min_size +
                      null_bitmap_size)

    print("\n--- Final Calculation ---")
    print("The minimum record size is the sum of the non-nullable data and the null bitmap overhead.")
    print(f"Total Size = Size('flightNumber') + Size('ticketCost') + Size('arrivalCity') + Size('departureCity') + Size(Null Bitmap)")
    print(f"Total Size = {flight_number_min_size} + {ticket_cost_min_size} + {arrival_city_min_size} + {departure_city_min_size} + {null_bitmap_size}")
    print(f"Total Minimum Storage Required: {total_min_size} bytes")

    # Return the final numeric answer for the platform.
    return total_min_size

# Execute the calculation and store the final answer.
final_answer = calculate_minimum_record_size()
print(f"\n<<<5>>>")
