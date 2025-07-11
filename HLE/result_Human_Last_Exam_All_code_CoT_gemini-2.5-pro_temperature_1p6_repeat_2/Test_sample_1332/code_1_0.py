import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """

    # --- Step 1: Calculate storage for the non-nullable PRIMARY KEY ---
    # The 'flightNumber' is an integer and a PRIMARY KEY, so it cannot be NULL.
    # A standard integer typically requires 4 bytes.
    flight_number_bytes = 4
    print(f"Storage for 'flightNumber' (integer, NOT NULL): {flight_number_bytes} bytes")

    # --- Step 2: Calculate storage for nullable fields in their smallest state ---
    # To find the minimum record size, we assume all nullable fields are NULL.
    # When a field is NULL, it doesn't store its data type's value (e.g., 4 bytes for an int).
    # Its data storage size is effectively 0 bytes.
    ticket_cost_bytes = 0
    arrival_city_bytes = 0
    departure_city_bytes = 0
    print(f"Storage for 'ticketCost' (when NULL): {ticket_cost_bytes} bytes")
    print(f"Storage for 'arrivalCity' (when NULL): {arrival_city_bytes} bytes")
    print(f"Storage for 'departureCity' (when NULL): {departure_city_bytes} bytes")
    
    # --- Step 3: Account for the Null Bitmap overhead ---
    # The database needs to track which of the 4 columns are NULL.
    # This is done with a null bitmap where each column gets one bit.
    # 4 columns require 4 bits. Since storage is byte-aligned, this
    # is rounded up to the nearest byte.
    total_columns = 4
    null_bitmap_bytes = math.ceil(total_columns / 8)
    print(f"Overhead for Null Bitmap ({total_columns} columns): {null_bitmap_bytes} byte(s)")

    # --- Step 4: Calculate the total minimum size ---
    total_min_bytes = (flight_number_bytes +
                       ticket_cost_bytes +
                       arrival_city_bytes +
                       departure_city_bytes +
                       null_bitmap_bytes)

    # --- Step 5: Print the final equation ---
    print("\n--- Final Calculation ---")
    print("The minimum storage space is the sum of the space for non-nullable data and the null bitmap.")
    print(f"Equation: {flight_number_bytes} (flightNumber) + {ticket_cost_bytes} (ticketCost) + {arrival_city_bytes} (arrivalCity) + {departure_city_bytes} (departureCity) + {null_bitmap_bytes} (null_bitmap) = {total_min_bytes} bytes")
    print(f"\nTotal Minimum Storage Required: {total_min_bytes} bytes")


calculate_min_record_size()
<<<5>>>