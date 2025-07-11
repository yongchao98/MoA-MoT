def calculate_minimum_storage():
    """
    Calculates the minimum storage space in bytes for a single record in the
    FlightRecords table.

    The calculation is based on the following assumptions for minimum size:
    - The primary key 'flightNumber' (integer) cannot be NULL.
    - All other fields ('ticketCost', 'arrivalCity', 'departureCity') are NULL.
    """
    
    # Standard size for an integer data type in bytes
    size_of_integer = 4
    
    # Size for a NULL value (data portion is 0 bytes)
    size_of_null = 0
    
    # Field sizes for the smallest possible record
    flightnumber_size = size_of_integer
    ticketcost_size = size_of_null
    arrivalcity_size = size_of_null
    departurecity_size = size_of_null
    
    # Calculate the total minimum size
    total_min_size = flightnumber_size + ticketcost_size + arrivalcity_size + departurecity_size
    
    print("Analysis for Minimum Record Size:")
    print(f" - 'flightNumber' (integer, PRIMARY KEY): Must have a value. Size = {flightnumber_size} bytes.")
    print(f" - 'ticketCost' (integer, nullable): Assumed to be NULL for minimum size. Size = {ticketcost_size} bytes.")
    print(f" - 'arrivalCity' (text, nullable): Assumed to be NULL for minimum size. Size = {arrivalcity_size} bytes.")
    print(f" - 'departureCity' (text, nullable): Assumed to be NULL for minimum size. Size = {departurecity_size} bytes.")
    print("\nMinimum required storage space is the sum of the individual field sizes.")
    print(f"Final Equation: {flightnumber_size} + {ticketcost_size} + {arrivalcity_size} + {departurecity_size} = {total_min_size} bytes")

calculate_minimum_storage()