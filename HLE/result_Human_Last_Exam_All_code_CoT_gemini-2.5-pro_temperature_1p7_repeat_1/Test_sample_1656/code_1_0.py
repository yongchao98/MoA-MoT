def find_permutation_cycles(p, n):
    """
    Finds the number of cycles in a permutation.
    The permutation is given as a dictionary mapping i to p(i).
    n is the size of the set {1, ..., n}.
    """
    visited = [False] * (n + 1)
    cycles = 0
    for i in range(1, n + 1):
        if not visited[i]:
            cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = p[j]
    return cycles

def solve_braid_index():
    """
    Solves for the braid index of the knot corresponding to the given grid diagram.
    """
    n = 7
    o_coords = [(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)]
    # x_coords are not needed for this calculation but are part of the problem description.
    # x_coords = [(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)]

    # The permutation pi_O maps column i to row j.
    pi_o = {i: j for i, j in o_coords}
    
    # Calculate the number of cycles in pi_O
    c_pi_o = find_permutation_cycles(pi_o, n)

    # The braid index can be calculated using the formula: n - c(pi_O)
    braid_index = n - c_pi_o
    
    print("The grid number is n = {}".format(n))
    print("The permutation pi_O, derived from the o-markings, is: {}".format(pi_o))
    print("The number of cycles in pi_O is c(pi_O) = {}".format(c_pi_o))
    print("The braid index is calculated as n - c(pi_O).")
    print("Braid Index = {} - {} = {}".format(n, c_pi_o, braid_index))

solve_braid_index()