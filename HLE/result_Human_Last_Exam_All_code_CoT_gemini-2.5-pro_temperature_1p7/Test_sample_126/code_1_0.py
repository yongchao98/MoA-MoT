import numpy as np

def print_gate(name, gate_matrix):
    """Helper function to print a gate matrix."""
    print(f"{name} =")
    for row in gate_matrix:
        # Format each element to be readable, handling complex numbers
        formatted_row = [f"{x:^4}" for x in row]
        print(" ".join(formatted_row))
    print()

def get_swap_variant(A=1, F=0, G=1, J=1, K=0, P=1):
    """
    Creates a 4x4 matrix based on the structure of Option D.
    Default values correspond to the standard SWAP gate.
    """
    # Ensure inputs are complex numbers for calculations
    A, F, G, J, K, P = map(complex, [A, F, G, J, K, P])
    
    # Structure from Option D
    gate = np.array([
        [A, 0, 0, 0],
        [0, F, G, 0],
        [0, J, K, 0],
        [0, 0, 0, P]
    ])
    return gate

# --- Demonstration ---

# 1. Show the general structure of Option D is the parent form
# We use strings for the variables to show the symbolic form
print("The general matrix form for Option D is:")
general_form = [
    ['A', '0', '0', '0'],
    ['0', 'F', 'G', '0'],
    ['0', 'J', 'K', '0'],
    ['0', '0', '0', 'P']
]
print_gate("U(A,F,G,J,K,P)", general_form)

# 2. Show that the problem's "generalized form" is a special case where F=0, K=0
print("Setting F=0 and K=0 gives the problem's 'generalized form':")
generalized_form = get_swap_variant(A='A', F=0, G='G', J='J', K=0, P='P')
print_gate("U(A,G,J,P)", [['A', 0, 0, 0], [0, 0, 'G', 0], [0, 'J', 0, 0], [0, 0, 0, 'P']])


# 3. Show that specific examples from the prompt can be created from this structure.

# Standard SWAP gate: A=1, P=1, F=0, K=0, G=1, J=1
swap_gate = get_swap_variant(A=1, F=0, G=1, J=1, K=0, P=1)
print_gate("SWAP Gate", swap_gate)

# iSWAP gate: A=1, P=1, F=0, K=0, G=i, J=i
# In Python, complex numbers are written with 'j'
iswap_gate = get_swap_variant(A=1, F=0, G=1j, J=1j, K=0, P=1)
print_gate("iSWAP Gate", iswap_gate)

# fSWAP gate: A=1, P=1, F=0, K=0, G=-1, J=-1
fswap_gate = get_swap_variant(A=1, F=0, G=-1, J=-1, K=0, P=1)
print_gate("fSWAP Gate", fswap_gate)
