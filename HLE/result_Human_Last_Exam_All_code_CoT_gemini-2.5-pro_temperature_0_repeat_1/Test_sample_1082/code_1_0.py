def calculate_sva():
    """
    Calculates the Sustainable Value Added (SVA) based on the provided details.

    The SVA is calculated as the revenue generated minus the costs associated
    with the eco-efficient (sustainable) benchmark.
    """
    # 1. Identify the revenue from the sale of the product.
    revenue = 50  # in dollars

    # 2. Identify the costs for the optimal/sustainable use of resources.
    cost_water = 10  # in dollars
    cost_energy = 15  # in dollars

    # 3. Calculate the total sustainable cost.
    total_sustainable_cost = cost_water + cost_energy

    # 4. Calculate the Sustainable Value Added.
    sva = revenue - total_sustainable_cost

    # 5. Print the final equation and the result.
    print(f"Sustainable Value Added = ${revenue} (Revenue) - (${cost_water} (Water Cost) + ${cost_energy} (Energy Cost)) = ${sva}")

calculate_sva()