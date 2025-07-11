def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided financial figures.

    The function interprets the problem's financial data as:
    - Revenue = $50
    - Material Cost = $10
    - Water Cost = $10
    - Energy Cost = $15

    Other quantitative data (250g, 100L, 40kWh, etc.) are treated as
    contextual information, as there is insufficient data to link them
    to a cost-per-unit calculation.

    SVA is calculated as Revenue - Total Costs.
    The final output prints the full equation.
    """
    revenue = 50
    material_cost = 10
    water_cost = 10
    energy_cost = 15

    # Calculate total cost and SVA
    total_cost = material_cost + water_cost + energy_cost
    sva = revenue - total_cost

    # Print the final equation with all components, as requested.
    print(f"The Sustainable Value Added (SVA) is calculated as Revenue - Total Costs.")
    print(f"Using the provided values:")
    print(f"SVA = ${revenue} - (${material_cost} + ${water_cost} + ${energy_cost})")
    print(f"SVA = ${revenue} - ${total_cost}")
    print(f"Final SVA = ${sva}")


calculate_sva()