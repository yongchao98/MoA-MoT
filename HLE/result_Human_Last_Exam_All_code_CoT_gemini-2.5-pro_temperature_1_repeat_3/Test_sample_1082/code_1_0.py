def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided financial data.
    """
    # Define the given financial values from the problem description
    selling_price = 50
    water_cost = 10
    energy_cost = 15

    # The Sustainable Value Added (SVA) is the value created after accounting for environmental costs.
    # It is calculated as the revenue (Value Added) minus the cost of environmental impacts.
    
    # Calculate the total environmental cost
    total_environmental_cost = water_cost + energy_cost
    
    # Calculate the final SVA
    sva = selling_price - total_environmental_cost
    
    # Print the explanation and the final equation with all the numbers, as requested.
    print("The formula for Sustainable Value Added (SVA) is: Revenue - Environmental Costs")
    print("\nCalculating the SVA for the product:")
    print(f"Selling Price = ${selling_price}")
    print(f"Water Cost = ${water_cost}")
    print(f"Energy Cost = ${energy_cost}")
    
    print("\nThe final equation is:")
    print(f"SVA = {selling_price} - ({water_cost} + {energy_cost})")
    
    print(f"\nResult:")
    print(f"SVA = {selling_price} - {total_environmental_cost} = ${sva}")

calculate_sustainable_value_added()