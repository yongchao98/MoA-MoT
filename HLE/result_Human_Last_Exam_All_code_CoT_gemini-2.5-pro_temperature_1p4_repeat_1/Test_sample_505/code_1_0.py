import sys

def solve_qubit_encoding():
    """
    Calculates the maximum number of logical qubits for a given surface code setup.

    The number of logical qubits 'k' on a single planar surface code patch is
    equal to the number of holes 'h' on that patch (k = h).

    For multiple, separate patches, the total number of logical qubits is the
    sum of the qubits on each patch. Therefore, the total number of logical
    qubits is equal to the total number of holes.
    """
    # Total number of holes specified in the problem.
    total_holes = 2
    
    # The problem specifies two separate patches. The distribution of the 2 holes
    # between these patches does not change the total number of logical qubits.
    # We can illustrate with any valid distribution, e.g., 1 hole per patch.
    holes_patch1 = 1
    holes_patch2 = 1
    
    # For each patch, the number of logical qubits equals its number of holes.
    logical_qubits_patch1 = holes_patch1
    logical_qubits_patch2 = holes_patch2
    
    # The total number of logical qubits is the sum from each patch.
    total_logical_qubits = logical_qubits_patch1 + logical_qubits_patch2
    
    print("The total number of logical qubits is the sum of logical qubits from each patch.")
    print("Let's demonstrate with a distribution of 1 hole per patch:")
    print(f"Logical qubits from Patch 1 (with {holes_patch1} hole): {logical_qubits_patch1}")
    print(f"Logical qubits from Patch 2 (with {holes_patch2} hole): {logical_qubits_patch2}")
    
    print("\nFinal Equation:")
    # The final output needs to show each number in the equation.
    print(f"{logical_qubits_patch1} + {logical_qubits_patch2} = {total_logical_qubits}")
    
    print(f"\nThus, at most {total_logical_qubits} logical qubits can be encoded.")

solve_qubit_encoding()

# The final answer as a raw value.
sys.stdout.write("<<<2>>>\n")