def solve_extinction_rate_factor():
    """
    This script solves for the multiplicative factor by which the extinction rate
    of a morphospecies is greater than that of an evolutionary species.
    It follows a step-by-step logical derivation based on the problem statement.
    """

    print("Let's define the rates of the underlying biological processes:")
    print("  - mu: the rate of true extinction of an evolutionary lineage.")
    print("  - lambda: the rate of true branching speciation of a lineage.")
    print("  - a: the rate of anagenetic change where one morphospecies is replaced by a new one in the same lineage.")
    print("-" * 50)

    # Step 1: Define the extinction rate for an Evolutionary Species (ES)
    print("Step 1: Extinction Rate of an Evolutionary Species (mu_es)")
    print("An evolutionary species goes extinct only when its lineage dies. This occurs at rate mu.")
    print("mu_es = 1 * mu")
    print("-" * 50)

    # Step 2: Define the extinction rate for a Morphospecies (MS)
    print("Step 2: Extinction Rate of a Morphospecies (mu_ms)")
    print("A morphospecies goes extinct from true extinction or pseudo-extinction.")
    print("The components of its extinction rate are:")
    print("  - True Extinction: The lineage dies out (rate: mu).")
    print("  - Pseudo-extinction via Anagenesis: The species is renamed (rate: a).")
    print("  - Pseudo-extinction via Bifurcation: The mother species is replaced by two daughters. This occurs in 50% of branching events (rate: 0.5 * lambda).")
    print("Therefore, the total extinction rate for a morphospecies is the sum of these rates:")
    print("mu_ms = mu + a + 0.5 * lambda")
    print("-" * 50)

    # Step 3: Apply the "same rates" assumption
    print("Step 3: Interpreting the 'Same Rates' Assumption")
    print("The problem states 'all the processes that affect them occur at the same rates'.")
    print("We can categorize extinction into two fundamental processes: 'true extinction' and 'pseudo-extinction'.")
    print("  - Rate of True Extinction = mu")
    print("  - Rate of Pseudo-extinction = a + 0.5 * lambda")
    print("We apply the assumption by equating the rates of these two processes:")
    print("Rate(True Extinction) = Rate(Pseudo-extinction)")
    print("mu = a + 0.5 * lambda")
    print("-" * 50)

    # Step 4: Calculate the final ratio
    print("Step 4: Calculating the Final Multiplicative Factor")
    print("We want to find the ratio mu_ms / mu_es.")
    print("Let's take the formula for mu_ms and group the pseudo-extinction terms:")
    print("mu_ms = mu + (a + 0.5 * lambda)")
    print("From our assumption in Step 3, we can substitute 'mu' for '(a + 0.5 * lambda)':")
    print("mu_ms = mu + mu")
    print("mu_ms = 2 * mu")
    print("\nNow we can write the final equation for the ratio:")
    print("Ratio = mu_ms / mu_es")
    # The final equation with each number outputted
    print("Ratio = (2 * mu) / (1 * mu)")
    final_factor = 2 / 1
    print(f"\nThe mu terms cancel out, leaving the final factor.")

solve_extinction_rate_factor()

<<<2>>>