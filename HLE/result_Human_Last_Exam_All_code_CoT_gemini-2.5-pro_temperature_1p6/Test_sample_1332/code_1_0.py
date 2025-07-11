import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a single record
    in the FlightRecords table.
    """

    # --- Step 1: Calculate size of non-nullable fields ---
    # The 'flightNumber' is an integer PRIMARY KEY, so it cannot be NULL.
    # A standard integer takes up 4 bytes.
    size_of_integer = 4
    flightnumber_size = size_of_integer
    print(f"Analyzing the storage requirements for the smallest possible record:")
    print(f"1. The 'flightNumber' (integer) is a PRIMARY KEY and cannot be null.")
    print(f"   - Storage for flightNumber: {flightnumber_size} bytes")

    # --- Step 2: Account for nullable fields ---
    # To have the smallest possible record, all nullable fields will be NULL.
    # A NULL value itself takes up 0 bytes in the data area.
    print(f"2. The 'ticketCost', 'arrivalCity', and 'departureCity' fields can be NULL.")
    print(f"   - To minimize size, these will be NULL, contributing 0 bytes of data.")

    # --- Step 3: Calculate the NULL bitmap overhead ---
    # Databases use a NULL bitmap in the row header to track null values.
    # This bitmap requires 1 bit for each column.
    num_columns = 4  # (flightNumber, ticketCost, arrivalCity, departureCity)
    # The size of the bitmap is the number of bits, rounded up to the nearest byte.
    bitmap_size_bytes = math.ceil(num_columns / 8)
    print(f"3. An overhead is required to track which of the {num_columns} columns are NULL.")
    print(f"   - This is typically a NULL bitmap requiring 1 bit per column, for a total of {num_columns} bits.")
    print(f"   - The storage for the bitmap is rounded up to the nearest byte: ceil({num_columns}/8) = {bitmap_size_bytes} byte.")

    # --- Step 4: Calculate the total minimum size ---
    total_min_size = flightnumber_size + bitmap_size_bytes
    print("\n-------------------------------------------------------------")
    print("The minimum storage space is the sum of the non-nullable data and the NULL bitmap.")
    print(f"Final Equation: {flightnumber_size} (for flightNumber) + {bitmap_size_bytes} (for NULL bitmap) = {total_min_size} bytes")
    print("-------------------------------------------------------------")
    
    # This print is for the final answer extraction.
    print(f"\n<<<{total_min_size}>>>")

calculate_min_record_size()