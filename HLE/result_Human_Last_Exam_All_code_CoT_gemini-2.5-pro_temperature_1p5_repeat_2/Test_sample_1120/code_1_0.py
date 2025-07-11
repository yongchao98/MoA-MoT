def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price based on the problem description.
    """
    # Step 1: Determine the equilibrium quantity supplied.
    # The producer's profit is maximized by selling all 10 units, as MC=0.
    Q_supply = 10
    
    # In equilibrium, quantity demanded equals quantity supplied.
    Q_equilibrium = 10

    # Step 2: Calculate the equilibrium price (P) using the inverse demand function.
    # The inverse demand function is 100*P = 400 + 3*Q^2 - Q^3/20.
    # We substitute Q = 10 to find the price that clears the market.
    numerator_for_price = 400 + 3 * (Q_equilibrium**2) - (Q_equilibrium**3) / 20
    P_equilibrium = numerator_for_price / 100

    print(f"The profit-maximizing quantity for the producer to supply is {Q_supply} units.")
    print(f"At market equilibrium, the quantity demanded must also be {Q_equilibrium} units.")
    print(f"The price that results in a demand of {Q_equilibrium} units is calculated as follows:")
    print(f"100*P = 400 + 3*({Q_equilibrium})^2 - ({Q_equilibrium})^3/20 = {numerator_for_price}")
    print(f"P = {numerator_for_price} / 100 = {P_equilibrium}")
    print(f"\nThus, the equilibrium price is {P_equilibrium}.")
    
    # Step 3: Calculate the excess demand at the equilibrium price.
    # By definition, at equilibrium, Q_demand equals Q_supply.
    Q_demand_at_equilibrium = Q_equilibrium
    excess_demand = Q_demand_at_equilibrium - Q_supply
    
    print("\nExcess demand is the difference between quantity demanded and quantity supplied at the equilibrium price.")
    print("The final calculation is:")
    # The prompt requires printing each number in the final equation.
    print(f"{int(Q_demand_at_equilibrium)} - {int(Q_supply)} = {int(excess_demand)}")

solve_excess_demand()
<<<0>>>