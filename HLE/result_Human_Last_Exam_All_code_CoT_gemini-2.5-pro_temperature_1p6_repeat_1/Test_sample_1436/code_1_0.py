import math

def calculate_symmetry_factor():
    """
    Calculates the symmetry factor for a Feynman diagram with two vertices
    connected by four propagators.
    """
    # Number of identical propagators connecting the two vertices
    n_propagators = 4
    
    # Number of identical vertices that can be interchanged
    n_vertices = 2

    # The symmetry factor S is the product of the factorials of the counts of identical components.
    # 1. Symmetry from permuting the identical propagators:
    propagator_perms = math.factorial(n_propagators)
    
    # 2. Symmetry from permuting the identical vertices:
    vertex_perms = math.factorial(n_vertices)
    
    # Total symmetry factor is the product of these individual factors
    symmetry_factor = propagator_perms * vertex_perms

    print("This script calculates the symmetry factor 'S' for a Feynman diagram.")
    print("The diagram consists of 2 identical vertices connected by 4 identical propagators.\n")
    print("The formula for the symmetry factor is:")
    print("S = (permutations of identical propagators) * (permutations of identical vertices)\n")

    print(f"Number of identical propagators = {n_propagators}")
    print(f"Number of permutations for propagators = {n_propagators}! = {propagator_perms}\n")
    
    print(f"Number of identical vertices = {n_vertices}")
    print(f"Number of permutations for vertices = {n_vertices}! = {vertex_perms}\n")

    print("The final calculation is:")
    print(f"S = {propagator_perms} * {vertex_perms}")
    print(f"S = {symmetry_factor}")

if __name__ == "__main__":
    calculate_symmetry_factor()