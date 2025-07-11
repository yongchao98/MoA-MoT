import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record in the specified table.

    The calculation is based on the following logic:
    1.  Non-nullable fields always occupy space. 'flightNumber' is an INTEGER and PRIMARY KEY,
        so it cannot be null. A standard integer requires 4 bytes.
    2.  Nullable fields ('ticketCost', 'arrivalCity', 'departureCity') are assumed to be NULL
        to achieve the minimum size. The data for these fields will occupy 0 bytes.
    3.  A null bitmap is used to track which columns are NULL. It needs 1 bit per column.
        For 4 columns, that's 4 bits. Storage is allocated in whole bytes, so we round up
        to the nearest byte. ceil(4 bits / 8 bits_per_byte) = 1 byte.
    4.  The total minimum size is the sum of the storage for non-nullable data and the
        null bitmap.
    """
    # Size of a standard integer in bytes
    integer_size = 4

    # Number of columns in the table
    num_columns = 4  # flightNumber, ticketCost, arrivalCity, departureCity

    # --- Calculation ---

    # 1. Storage for non-nullable data
    # Only 'flightNumber' is non-nullable.
    non_nullable_data_size = integer_size

    # 2. Storage for the null bitmap
    # The bitmap needs 1 bit per column. The size is rounded up to the nearest full byte.
    bits_per_byte = 8
    null_bitmap_size = math.ceil(num_columns / bits_per_byte)

    # 3. Total minimum size
    total_min_size = non_nullable_data_size + null_bitmap_size

    # --- Output ---
    print("Calculating the minimum storage space for a single record:")
    print(f"- The 'flightNumber' (integer) is non-nullable and requires {non_nullable_data_size} bytes.")
    print("- All other fields are assumed to be NULL, so their data size is 0 bytes.")
    print(f"- A null bitmap is required to track nulls for all {num_columns} columns, which takes {null_bitmap_size} byte.")
    print("\nFinal equation for minimum storage space (in bytes):")
    print(f"size('flightNumber') + size(null_bitmap)")
    print(f"{non_nullable_data_size} + {null_bitmap_size} = {total_min_size}")


if __name__ == '__main__':
    calculate_min_record_size()