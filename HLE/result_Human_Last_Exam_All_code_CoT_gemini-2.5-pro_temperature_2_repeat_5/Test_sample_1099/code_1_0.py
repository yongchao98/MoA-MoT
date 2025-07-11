import math

def calculate_simulation_resources():
    """
    Calculates and explains the resources needed to simulate a singlet state
    with an LHV model.
    """
    
    # Part 1: Simulation cost using a PR-Box (no communication)

    # Define the maximal CHSH inequality values for the different models
    chsh_local = 2.0
    chsh_pr = 4.0
    chsh_quantum = 2 * math.sqrt(2)
    
    # To find the minimal fraction 'p' of a PR-Box needed, we solve the equation:
    # p * chsh_pr + (1-p) * chsh_local = chsh_quantum
    # p * 4 + (1-p) * 2 = 2*sqrt(2)
    # 4p + 2 - 2p = 2*sqrt(2)
    # 2p = 2*sqrt(2) - 2
    # p = sqrt(2) - 1
    
    p = math.sqrt(2) - 1

    print("To classically simulate the correlations of a singlet quantum state, a Local Hidden Variable (LHV) model must be supplemented with additional resources.")
    print("\n--- Resource Costs ---")
    
    print("\n1. Cost in terms of PR-Box non-locality:")
    print("This is the minimal fraction 'p' of a PR-Box which, when mixed with a local resource, can reproduce the quantum correlations.")
    print("The governing equation relates the CHSH Bell inequality values:")
    print("p * (CHSH_PR_Box) + (1-p) * (CHSH_Local) = CHSH_Quantum")
    print("\nSubstituting the values into the equation:")
    print(f"p * {chsh_pr} + (1-p) * {chsh_local} = {chsh_quantum:.8f}")
    
    print("\nSolving for p gives:")
    print(f"p = (CHSH_Quantum - CHSH_Local) / (CHSH_PR_Box - CHSH_Local)")
    print(f"p = ({chsh_quantum:.8f} - {chsh_local}) / ({chsh_pr} - {chsh_local})")
    print(f"p = {p:.8f}")
    print("\nThis cost is exactly sqrt(2) - 1, which represents the 'amount' of non-locality required.")

    print("\n-------------------------------------")

    # Part 2: Simulation cost using classical communication (no PR-Box)
    communication_cost_bits = 1.0

    print("\n2. Cost in terms of classical communication:")
    print("It has been shown (Toner and Bacon, 2003) that the same correlations can be simulated using an LHV model supplemented with, on average, a specific amount of communication.")
    print(f"Average communication cost: {communication_cost_bits} bit")

if __name__ == '__main__':
    calculate_simulation_resources()
