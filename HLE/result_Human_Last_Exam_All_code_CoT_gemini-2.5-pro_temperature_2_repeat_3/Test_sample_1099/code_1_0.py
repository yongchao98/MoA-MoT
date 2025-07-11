import math

def calculate_simulation_resources():
    """
    Calculates the minimal resources to simulate singlet state correlations.

    The problem is to find the minimal amount of non-local resources (PR-Box)
    and communication needed for a Local Hidden Variable (LHV) model to simulate
    the correlations from any POVM measurements on a singlet quantum state.

    The simulation cost is determined by the most non-local correlations,
    which are those that maximally violate the CHSH inequality.

    We use a simulation model that mixes local resources and PR-Box resources.
    Let 'c' be the fraction of PR-Box resource used on average.
    The maximum CHSH value of the simulation is then a weighted average of the
    maximum local CHSH value and the maximum PR-Box CHSH value.

    The governing equation is:
    S_quantum = (1 - c) * S_local + c * S_pr_box
    """

    # Maximum CHSH value for a local theory (Bell's limit)
    S_local = 2.0

    # Maximum CHSH value for a singlet quantum state (Tsirelson's bound)
    S_quantum = 2 * math.sqrt(2)

    # Maximum CHSH value for a non-signaling PR-Box
    S_pr_box = 4.0

    print("Step 1: Define the CHSH inequality bounds for different models.")
    print(f"  - Maximum for Local Hidden Variable models (S_local): {S_local}")
    print(f"  - Maximum for a Quantum Singlet State (S_quantum): {S_quantum:.4f}")
    print(f"  - Maximum for a PR-Box (S_pr_box): {S_pr_box}\n")

    print("Step 2: Set up the simulation equation.")
    print("  We need to find the minimal fraction 'c' of a PR-Box such that:")
    print("  S_quantum = (1 - c) * S_local + c * S_pr_box\n")

    print("Step 3: Substitute the known values into the equation.")
    # The equation with numbers: S_quantum = (1-c)*2 + c*4
    print("  The equation to solve for 'c' is:")
    print(f"  {S_quantum:.4f} = (1 - c) * {S_local} + c * {S_pr_box}")
    print(f"  {S_quantum:.4f} = {S_local} - {S_local}*c + {S_pr_box}*c")
    print(f"  {S_quantum:.4f} - {S_local} = ({S_pr_box} - {S_local}) * c\n")


    # Solving for c:
    # S_quantum = 2 + 2*c
    # c = (S_quantum - 2) / 2
    c = (S_quantum - S_local) / (S_pr_box - S_local)

    # The communication cost in this type of simulation is 0.
    communication_cost = 0

    print("Step 4: Solve for 'c' and state the final resources.")
    print(f"  c = ({S_quantum:.4f} - {S_local}) / ({S_pr_box} - {S_local})")
    print(f"  c = {c:.4f}\n")

    print("--- Final Answer ---")
    print("The minimal resources required on average are:")
    print(f"  - PR-Box Fraction (c): {c:.4f} (which is sqrt(2) - 1)")
    print(f"  - Classical Communication (bits): {communication_cost}")


if __name__ == '__main__':
    calculate_simulation_resources()
