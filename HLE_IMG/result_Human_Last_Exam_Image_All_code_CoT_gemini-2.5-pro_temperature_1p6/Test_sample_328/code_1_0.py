import collections

def solve_puzzle():
    """
    Solves the nonlinear wave equation plot identification puzzle.
    
    The logic behind the assignments is based on a detailed analysis of the equation's
    parameters (b, c, d) and their visual impact on the solution plots.
    - b: Controls symmetry (b=0 is symmetric, b>0 red-asymmetric, b<0 blue-asymmetric).
    - c: Controls stability/complexity (c=1 chaotic, c=-1 simple/regular, c=0 neutral).
    - d: Controls drift (d>0 red-drift, d<0 blue-drift, d=0 no drift).

    Each plot is matched to a parameter combination where one parameter is unique
    and the other two are equal, based on these visual cues.
    """
    
    # The solution is derived from a careful case-by-case analysis.
    # The assignments are stored in an ordered dictionary to preserve plot order.
    assignments = collections.OrderedDict()
    
    # Plot #: Code (Unique Param=Unique Val, Common Param=Common Val) # Visual Rationale
    assignments[1] = 'C'  # (c=1, b=d=0)   # Chaotic, symmetric, neutral-drift
    assignments[2] = 'b'  # (b=-1, c=d=1)  # Chaotic, blue-asym, red-drift
    assignments[3] = 'd'  # (d=-1, b=c=0)  # Symmetric, neutral-stab, blue-drift
    assignments[4] = 'b'  # (b=-1, c=d=0)  # Blue-asym, neutral-stab, no-drift
    assignments[5] = 'B'  # (b=1, c=d=-1)  # Red-asym, simple, blue-drift
    assignments[6] = 'D'  # (d=1, b=c=-1)  # Blue-asym, simple, red-drift
    assignments[7] = 'd'  # (d=-1, b=c=1)  # Red-asym, chaotic, blue-drift
    assignments[8] = 'Z'  # (c=0, b=d=1)   # Red-asym, neutral-stab, red-drift
    assignments[9] = 'c'  # (c=-1, b=d=1)  # Red-asym, simple, red-drift
    assignments[10] = 'c' # (c=-1, b=d=0)  # Symmetric, simple, no-drift
    assignments[11] = '0' # (d=0, b=c=1)   # Red-asym, chaotic, no-drift
    assignments[12] = 'z' # (b=0, c=d=-1)  # Symmetric, simple, blue-drift
    assignments[13] = 'Z' # (c=0, b=d=-1)  # Blue-asym, neutral-stab, blue-drift
    assignments[14] = '0' # (d=0, b=c=-1)  # Blue-asym, simple, no-drift
    assignments[15] = 'D' # (d=1, b=c=0)   # Symmetric, neutral-stab, red-drift
    assignments[16] = 'B' # (b=1, c=d=0)   # Red-asym, neutral-stab, no-drift
    
    # Construct the final answer string
    result_string = "".join(assignments.values())
    
    print(result_string)

solve_puzzle()