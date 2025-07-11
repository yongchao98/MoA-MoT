import math

def solve_optimal_mutation_rate():
    """
    This script models the trade-off that determines the optimal genomic mutation rate (U).
    It balances the fitness cost of deleterious mutations (mutational load)
    against the metabolic cost of maintaining high-fidelity DNA repair.

    - k1: A coefficient representing the fitness cost per mutation (mutational load).
    - k2: A coefficient representing the metabolic cost of DNA repair fidelity.
          A lower U implies higher fidelity and thus a higher cost.
    
    The optimal mutation rate U_optimal is the point where these costs are balanced,
    maximizing the overall fitness of the organism.
    """
    
    # Define the cost coefficients for the model
    k1 = 0.5  # Cost of mutational load
    k2 = 0.0005 # Cost of fidelity
    
    # Calculate the optimal mutation rate (U) that maximizes fitness
    # by solving for U in d(Fitness)/dU = 0, which yields U = sqrt(k2 / k1)
    U_optimal = math.sqrt(k2 / k1)
    
    print("Step 1: Define the cost coefficients.")
    print(f"Cost of mutational load (k1): {k1}")
    print(f"Cost of fidelity (k2): {k2}")
    print("\nStep 2: Calculate the optimal rate using the derived equation.")
    print("Equation: U_optimal = sqrt(k2 / k1)")
    # The instruction requested printing each number in the final equation.
    print(f"Calculation: U_optimal = sqrt({k2} / {k1})")
    
    print(f"\nResult: The calculated optimal genomic mutation rate (U_optimal) is approximately {U_optimal:.4f}.")
    print("\nThis model demonstrates how natural selection for fitness optimality (balancing opposing costs) can lead to a specific, constant genomic mutation rate.")

solve_optimal_mutation_rate()