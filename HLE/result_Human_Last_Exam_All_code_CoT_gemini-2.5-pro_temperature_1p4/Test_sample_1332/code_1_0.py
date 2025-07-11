def calculate_minimum_storage():
    """
    Calculates the minimum storage space in bytes for a record in the FlightRecords table.
    """
    # Step 1: Define the size of the data types in bytes.
    # A standard SQL INTEGER typically occupies 4 bytes.
    # A NULL value for any data type occupies 0 bytes of data storage.
    size_of_integer = 4
    size_of_null = 0

    # Step 2: Identify mandatory and nullable fields.
    # 'flightNumber' is a PRIMARY KEY, so it cannot be NULL. Its size is mandatory.
    flight_number_size = size_of_integer
    
    # 'ticketCost', 'arrivalCity', and 'departureCity' can be NULL.
    # To find the minimum storage, we assume they are all NULL.
    ticket_cost_size = size_of_null
    arrival_city_size = size_of_null
    departure_city_size = size_of_null
    
    # Step 3: Calculate the total minimum size by summing the sizes of each field.
    total_min_storage = flight_number_size + ticket_cost_size + arrival_city_size + departure_city_size
    
    print("This script calculates the minimum storage space based on non-nullable fields.")
    print("\n--- Equation Breakdown ---")
    print(f"Size of 'flightNumber' (integer, non-nullable): {flight_number_size} bytes")
    print(f"Size of 'ticketCost' (integer, nullable): {ticket_cost_size} bytes")
    print(f"Size of 'arrivalCity' (text, nullable): {arrival_city_size} bytes")
    print(f"Size of 'departureCity' (text, nullable): {departure_city_size} bytes")
    
    print("\n--- Final Calculation ---")
    print(f"{flight_number_size} + {ticket_cost_size} + {arrival_city_size} + {departure_city_size} = {total_min_storage}")
    
    print(f"\nTotal Minimum Storage Required: {total_min_storage} bytes")

calculate_minimum_storage()