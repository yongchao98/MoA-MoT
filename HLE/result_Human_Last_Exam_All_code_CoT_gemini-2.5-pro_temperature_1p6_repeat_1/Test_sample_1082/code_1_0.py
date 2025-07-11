def calculate_sustainable_value_added():
    """
    Calculates the Sustainable Value Added (SVA) based on provided financial figures.

    SVA is determined by subtracting the total environmental costs from the revenue generated.
    """
    
    # Revenue from selling the final product
    revenue = 50
    
    # Environmental costs for the production process
    water_cost = 10
    energy_cost = 15
    
    # Calculate the total environmental cost
    total_environmental_cost = water_cost + energy_cost
    
    # Calculate the Sustainable Value Added (SVA)
    sva = revenue - total_environmental_cost
    
    # Print the final equation as requested
    print(f"The Sustainable Value Added (SVA) is calculated as follows:")
    print(f"Revenue - (Water Cost + Energy Cost) = SVA")
    print(f"{revenue} - ({water_cost} + {energy_cost}) = {sva}")

# Execute the function
calculate_sustainable_value_added()