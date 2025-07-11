import math

def calculate_min_record_size():
    """
    Calculates the minimum storage space in bytes for a record based on a given table structure.

    This calculation considers:
    1. The fixed size of non-nullable columns.
    2. The space required for a null bitmap to track nullable columns.
    3. Assumes all nullable columns are set to NULL for the smallest possible record.
    4. Excludes database engine-specific overhead (e.g., row headers) for a general answer.
    """

    # Define the table structure and standard data type sizes in bytes.
    # A PRIMARY KEY is always non-nullable.
    columns = {
        'flightNumber': {'type': 'integer', 'nullable': False, 'size': 4},
        'ticketCost': {'type': 'integer', 'nullable': True, 'size': 4},
        'arrivalCity': {'type': 'text', 'nullable': True, 'size': 'variable'},
        'departureCity': {'type': 'text', 'nullable': True, 'size': 'variable'}
    }

    # --- Calculation ---

    # 1. Calculate space for fields that cannot be NULL.
    non_nullable_data_size = 0
    non_nullable_field_name = ""
    for name, properties in columns.items():
        if not properties['nullable']:
            # This field must contain data.
            non_nullable_data_size += properties['size']
            non_nullable_field_name = name

    # 2. Count the number of nullable fields to determine bitmap size.
    nullable_field_count = 0
    nullable_field_names = []
    for name, properties in columns.items():
        if properties['nullable']:
            nullable_field_count += 1
            nullable_field_names.append(name)
    
    # 3. Calculate null bitmap size. It needs 1 bit per nullable column, rounded up to the nearest byte.
    # Using integer arithmetic: (nullable_field_count + 7) // 8
    null_bitmap_size = (nullable_field_count + 7) // 8

    # 4. Sum the sizes for the final result.
    minimum_record_size = non_nullable_data_size + null_bitmap_size

    # --- Output ---
    
    print("To find the minimum storage space, we calculate the size of required data plus the overhead for tracking NULLs.")
    
    print("\n1. Size of Non-Nullable Data:")
    print(f"   - The '{non_nullable_field_name}' field is a PRIMARY KEY and cannot be NULL.")
    print(f"   - An 'integer' requires {columns[non_nullable_field_name]['size']} bytes.")
    print(f"   - Total size for required data = {non_nullable_data_size} bytes.")
    
    print("\n2. Size of Null-Tracking Overhead (Null Bitmap):")
    print(f"   - The table has {nullable_field_count} nullable fields: {', '.join(f'`{name}`' for name in nullable_field_names)}.")
    print(f"   - To track these, a null bitmap needs 1 bit per field, for a total of {nullable_field_count} bits.")
    print(f"   - Storage is byte-aligned, so we round up: ceil({nullable_field_count}/8) = {null_bitmap_size} byte(s).")

    print("\n3. Total Minimum Record Size:")
    print("   - This is the sum of the required data size and the null bitmap size.")
    print(f"   - Final Equation: {non_nullable_data_size} + {null_bitmap_size} = {minimum_record_size}")

    return minimum_record_size

if __name__ == "__main__":
    final_answer = calculate_min_record_size()
    print(f"\n<<<{final_answer}>>>")