import math
import itertools

def find_max_equilibria():
    """
    This function demonstrates that 2^N equilibria can exist for the given
    Lotka-Volterra system by considering a special case.
    
    The maximum number of equilibria is achieved when parameters are chosen
    such that the complex interaction term vanishes. This occurs when all
    A_i parameters are equal. The system then decouples into N independent
    logistic equations.
    """
    
    # The user can set the number of species (N) for the demonstration.
    N = 4

    # We choose parameters that satisfy the problem constraints (>=1).
    # For this demonstration, the specific values of R_i and Gamma don't affect
    # the location of the equilibria, only the dynamics. We need the K_i values.
    # Let's define the carrying capacities K_i. We make them distinct for clarity.
    K = [val + 2 for val in range(N)] 

    print(f"Demonstrating the number of equilibria for N = {N} species.")
    print("This maximum is achieved when all A_i parameters are equal.")
    print("In this case, each species' equilibrium abundance can be either 0 or its carrying capacity K_i.")
    print(f"Using carrying capacities K = {K}, the {2**N} possible equilibria are:\n")

    # For each species i, its equilibrium value can be 0 or K_i.
    # We generate all combinations of these values.
    equilibrium_options = [(0, k) for k in K]
    
    # itertools.product generates the Cartesian product, which corresponds to all
    # possible equilibrium points.
    equilibria_list = itertools.product(*equilibrium_options)
    
    count = 0
    for eq in equilibria_list:
        count += 1
        print(f"Equilibrium {count}: {eq}")

    print(f"\nTotal number of equilibria found: {count}")
    
    # The final equation for the maximum number of equilibria is 2^N.
    # We print the components of this equation as requested.
    base = 2
    exponent = N
    result = count # which is base**exponent
    print(f"\nThe final equation for the maximum number of equilibria is: {base}^{exponent} = {result}")

if __name__ == '__main__':
    find_max_equilibria()