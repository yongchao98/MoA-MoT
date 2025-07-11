import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the specified table.

    This calculation considers the size of the data itself plus the necessary
    overhead for tracking NULL values (a null bitmap). It ignores larger,
    database-system-specific overhead like full record headers or alignment padding.
    """

    # --- Step 1: Define sizes for data types ---
    # A standard integer requires 4 bytes. As a PRIMARY KEY, flightNumber cannot be NULL.
    flight_number_size = 4

    # For minimum record size, all nullable fields are assumed to be NULL.
    # A NULL value takes up 0 bytes for the data itself.
    ticket_cost_size_min = 0
    arrival_city_size_min = 0
    departure_city_size_min = 0

    # --- Step 2: Calculate overhead for tracking NULLs ---
    # A null bitmap is needed to track which columns are NULL.
    # It requires 1 bit per column. The size is rounded up to the nearest byte.
    num_columns = 4  # flightNumber, ticketCost, arrivalCity, departureCity
    null_bitmap_size = math.ceil(num_columns / 8.0)

    # --- Step 3: Sum the sizes to find the total minimum storage ---
    total_min_size = (flight_number_size +
                      ticket_cost_size_min +
                      arrival_city_size_min +
                      departure_city_size_min +
                      null_bitmap_size)

    # --- Step 4: Print the detailed explanation and final equation ---
    print("To find the minimum storage space, we calculate the size for a record where all nullable fields are NULL.")
    print("-" * 50)
    print(f"1. Size of Non-Nullable Data:")
    print(f"   - flightNumber (integer, NOT NULL): {flight_number_size} bytes")
    print("\n2. Size of Nullable Data (assumed NULL for minimum size):")
    print(f"   - ticketCost (integer, NULL): {ticket_cost_size_min} bytes")
    print(f"   - arrivalCity (text, NULL): {arrival_city_size_min} bytes")
    print(f"   - departureCity (text, NULL): {departure_city_size_min} bytes")
    print("\n3. Overhead for Null Bitmap:")
    print(f"   A bitmap is needed to track {num_columns} columns, requiring CEILING({num_columns}/8) bytes.")
    print(f"   - Null Bitmap Size: {null_bitmap_size} byte")
    print("-" * 50)
    print("Final Calculation:")
    print(f"Total Size = Size(flightNumber) + Size(ticketCost) + Size(arrivalCity) + Size(departureCity) + Size(Bitmap)")
    print(f"Total Size = {flight_number_size} + {ticket_cost_size_min} + {arrival_city_size_min} + {departure_city_size_min} + {null_bitmap_size} = {total_min_size} bytes")
    
    # Final answer in the required format
    print(f"\n<<<{total_min_size}>>>")

calculate_min_record_size()