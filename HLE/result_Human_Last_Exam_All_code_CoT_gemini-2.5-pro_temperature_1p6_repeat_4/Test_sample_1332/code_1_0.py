def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """

    # Size of a non-nullable integer (PRIMARY KEY) in bytes.
    flight_number_size = 4
    
    # Size of nullable fields when they are NULL is 0 bytes for the data itself.
    ticket_cost_size = 0
    arrival_city_size = 0
    departure_city_size = 0
    
    # Overhead to track NULL values.
    # A NULL bitmap is used, requiring 1 bit per column. With 4 columns, we need 4 bits.
    # Storage is allocated in whole bytes, so this rounds up to 1 byte.
    null_bitmap_size = 1
    
    # Total minimum size is the sum of the fixed data size and the NULL bitmap.
    total_size = flight_number_size + ticket_cost_size + arrival_city_size + departure_city_size + null_bitmap_size

    print("To determine the minimum storage space, we assume all optional fields are NULL.")
    print("-" * 30)
    print(f"flightNumber (integer, NOT NULL): {flight_number_size} bytes")
    print(f"ticketCost (set to NULL): {ticket_cost_size} bytes")
    print(f"arrivalCity (set to NULL): {arrival_city_size} bytes")
    print(f"departureCity (set to NULL): {departure_city_size} bytes")
    print(f"NULL Bitmap (to track NULLs): {null_bitmap_size} byte")
    print("-" * 30)
    print("Final Calculation:")
    print(f"Total = {flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} + {null_bitmap_size} = {total_size} bytes")

calculate_min_record_size()