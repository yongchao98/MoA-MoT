import sys

def solve_diversification_puzzle():
    """
    This script solves the diversification rate puzzle by defining and relating
    the extinction rates for evolutionary species and morphospecies.
    """

    # Step 1: Define the components of the morphospecies extinction rate (mu_m).
    # Based on the context, a morphospecies is considered "extinct" if one of the following occurs:
    # 1. The entire lineage goes extinct (true extinction, at rate mu_e).
    # 2. The lineage bifurcates, replacing the mother species with two daughter species.
    #    This happens 50% of the time a true speciation event occurs (rate = 0.5 * lambda_e).
    # 3. The paleontologist decides the morphology has changed enough to name a new species
    #    in the same lineage (anagenesis, at rate sigma), making the old one extinct.
    
    print("Step 1: Define the total extinction rate for a morphospecies (mu_m).")
    print("The total rate is the sum of rates from all events that cause a morphospecies to disappear.")
    print("mu_m = (rate of true extinction) + (rate of extinction from bifurcation) + (rate of extinction from anagenesis)")
    print("mu_m = mu_e + (0.5 * lambda_e) + sigma")
    print("-" * 30)

    # Step 2: Interpret the key assumption to relate the unknown rates.
    # The problem states: "Assume that for both evolutionary species and morphospecies, all the processes that affect them occur at the same rates."
    # We can interpret this as the rate of extinction due to biological reality (mu_e) being equal to the
    # rate of extinction due to taxonomic artifacts (pseudo-extinction).
    # The rate of pseudo-extinction is the sum of extinctions from bifurcation and anagenesis.
    
    print("Step 2: Apply the key assumption from the problem statement.")
    print("Assumption: Rate of true extinction = Rate of pseudo-extinction")
    print("mu_e = (0.5 * lambda_e) + sigma")
    print("-" * 30)

    # Step 3: Substitute the assumption into the equation for mu_m to find the final relationship.
    # We can replace the pseudo-extinction part of the mu_m formula with mu_e.
    
    print("Step 3: Calculate the relationship between mu_m and mu_e.")
    print("From Step 1, we have: mu_m = mu_e + (0.5 * lambda_e + sigma)")
    print("From Step 2, we substitute (0.5 * lambda_e + sigma) with mu_e:")
    print("mu_m = mu_e + mu_e")
    print("This simplifies to the final equation:")
    # The final equation shows that mu_m is twice mu_e.
    final_equation_lhs = "mu_m"
    factor = 2
    final_equation_rhs = "mu_e"
    print(f"{final_equation_lhs} = {factor} * {final_equation_rhs}")
    print("-" * 30)
    
    # Step 4: State the final answer.
    # The question asks how greater the extinction rate for a morphospecies is compared to an evolutionary species.
    # This is the multiplicative factor we just found.
    
    print("Step 4: Determine the multiplicative factor.")
    print(f"The extinction rate for a morphospecies (mu_m) is {factor} times greater than the extinction rate for an evolutionary species (mu_e).")
    
    # Required final output format
    sys.stdout.write(f"<<<{factor}>>>")

solve_diversification_puzzle()