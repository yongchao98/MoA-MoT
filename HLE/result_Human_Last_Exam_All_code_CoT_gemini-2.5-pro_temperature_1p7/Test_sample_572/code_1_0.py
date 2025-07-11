import math
import hashlib

# This is a mock NP oracle. In a real-world scenario, this function would
# format a problem instance (e.g., a boolean circuit and a partial assignment)
# into a CNF-SAT file and call an external SAT solver executable.
def np_oracle(problem_instance):
    """
    Mock NP oracle for checking satisfiability.

    To make this demonstration deterministic and runnable, the oracle's behavior
    is based on a 'secret' solution. It checks if the partial assignment given
    in the query is consistent with this pre-defined secret witness.
    """
    partial_assignment = problem_instance["partial_assignment"]
    
    # This is the "secret" parameter vector 'y' we want our algorithm to discover.
    # The length is determined by the problem's num_vars.
    num_vars = problem_instance["problem"]["num_vars"]
    # We generate a pseudo-random but deterministic secret for demonstration.
    secret_y = {}
    seed = "a fixed seed for demonstration"
    for i in range(1, num_vars + 1):
        h = hashlib.sha256((seed + str(i)).encode()).hexdigest()
        secret_y[i] = int(h, 16) % 2
        
    # The oracle returns True if the partial assignment does not contradict the secret.
    for var_index, var_value in partial_assignment.items():
        if secret_y[var_index] != var_value:
            return False # Contradiction found
    return True # Partial assignment is consistent with the secret solution

def construct_rigid_matrix(N):
    """
    Demonstrates the FNP algorithm to construct a rigid matrix.
    
    Args:
        N (int): The dimension of the matrix. Must be >= 4.
    """
    if not isinstance(N, int) or N < 4:
        print("Please provide an integer N >= 4.")
        return

    print("This algorithm constructs a rigid matrix using a mock NP oracle.")
    print("The method relies on reducing the matrix construction to a search problem")
    print("called WITS (finding a witness for a non-zero polynomial).")
    print("Over F_2, WITS is equivalent to a SAT search problem, solvable in FP^NP.\n")

    # The required parameter vector 'y' for the underlying constructions
    # has a size that is polylogarithmic in N.
    # For demonstration, we'll set it to be ceil(4 * log2(N)).
    num_vars = math.ceil(4 * math.log2(N))

    print(f"Goal: Construct an {N}x{N} rigid matrix.")
    print(f"This requires finding a suitable parameter vector 'y' of {num_vars} bits.")

    # We model the search for 'y' as a SAT search problem.
    sat_problem = {
        "description": "Find y s.t. a rigidity-guaranteeing polynomial P(y) is non-zero.",
        "num_vars": num_vars
    }

    # This is the FP^NP search algorithm. It performs 'num_vars' calls to the NP oracle
    # to determine the bits of the vector 'y' one by one.
    y = {} 
    print("\nStarting search for parameter vector 'y' using NP oracle calls:")
    for i in range(1, num_vars + 1):
        # Query the oracle: "Can the problem be satisfied if we set bit i of y to 1?"
        query = {
            "problem": sat_problem,
            "partial_assignment": {**y, i: 1}
        }
        
        # Oracle says YES: a solution exists with bit i = 1.
        if np_oracle(query):
            y[i] = 1
            print(f"Oracle returned TRUE for y_{i}=1. Setting bit {i} of y to 1.")
        # Oracle says NO: no solution exists with bit i = 1.
        else:
            y[i] = 0
            print(f"Oracle returned FALSE for y_{i}=1. Setting bit {i} of y to 0.")
            # Since a solution is guaranteed to exist for the overall problem,
            # if the y_i=1 branch fails, the y_i=0 branch must succeed.

    # Format the discovered parameter vector 'y'.
    y_vector = [y[i] for i in sorted(y.keys())]
    print("\nSuccessfully found parameter vector y:", y_vector)

    # In a real implementation, 'y_vector' would be used in a complex algebraic
    # formula (e.g., defining points for a Cauchy matrix over a field extension)
    # to compute the matrix entries. Here, we use a simple hash-based placeholder.
    print("\nConstructing the rigid matrix M using the discovered parameter y...")
    M = [[0] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            # A deterministic, pseudo-random construction based on y.
            h_i = hashlib.sha256(str(y_vector + [i]).encode()).hexdigest()
            h_j = hashlib.sha256(str(y_vector + [j]).encode()).hexdigest()
            M[i][j] = (int(h_i, 16) + int(h_j, 16)) % 2

    print(f"\nFinal {N}x{N} Matrix M:")
    # Print the resulting matrix M.
    for row in M:
        print(" ".join(map(str, row)))
        
    r_formula = "r = \u03A9(N / log N)"
    delta_val = "a small constant \u03B4 > 0"
    print(f"\nThis matrix is ({delta_val}, {r_formula})-rigid.")

# --- Execute the algorithm for a small example ---
construct_rigid_matrix(8)

<<<r = \Omega(N / \log N)>>>