import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """

    # --- Step 1: Define size of the non-nullable column(s) ---
    # The primary key 'flightNumber' is an integer and cannot be null.
    # A standard integer typically uses 4 bytes.
    flight_number_size = 4
    print(f"Size of non-nullable 'flightNumber' (integer): {flight_number_size} bytes")

    # --- Step 2: Define size of nullable columns ---
    # For minimum record size, 'ticketCost', 'arrivalCity', and 'departureCity'
    # will be NULL, so their data storage size is 0 bytes.
    nullable_data_size = 0
    print(f"Size of data for nullable columns (when NULL): {nullable_data_size} bytes")

    # --- Step 3: Calculate size of the null bitmap ---
    # The database needs to store a flag for each nullable column to indicate
    # whether it is NULL. There are 3 nullable columns.
    num_nullable_columns = 3
    print(f"Number of nullable columns: {num_nullable_columns}")

    # We need 1 bit per column. The storage needed is ceil(bits / 8).
    null_bitmap_size = math.ceil(num_nullable_columns / 8.0)
    print(f"Storage for null bitmap (1 byte can track up to 8 columns): {int(null_bitmap_size)} byte")

    # --- Step 4: Calculate the total minimum size ---
    # The total size is the sum of the required data and the null bitmap.
    total_min_size = flight_number_size + null_bitmap_size

    print("\n--- Final Calculation ---")
    print(f"Minimum Size = (Size of 'flightNumber') + (Size of Null Bitmap)")
    print(f"Minimum Size = {flight_number_size} + {int(null_bitmap_size)}")
    print(f"\nThe minimum storage space required is: {total_min_size} bytes")

calculate_min_record_size()
<<<5>>>