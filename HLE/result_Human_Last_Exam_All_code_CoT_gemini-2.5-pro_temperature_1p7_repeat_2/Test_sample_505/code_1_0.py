import math

def solve_logical_qubits():
    """
    Calculates the maximum number of logical qubits for the given surface code configuration.
    
    The problem specifies two patches of surface code with a total of two holes.
    """
    
    # Parameters from the problem statement
    num_patches = 2
    total_holes = 2
    
    # The number of logical qubits 'k' in a surface code is given by the formula: k = 2g + b - 1,
    # where 'g' is the genus and 'b' is the number of boundaries.
    
    # For a single planar "patch", the genus 'g' is 0. The number of boundaries 'b' is the
    # number of holes 'h' plus 1 for the outer boundary (b = h + 1).
    # So, for a single patch: k_patch = 2*0 + (h + 1) - 1 = h.
    # This means the number of logical qubits in a patch is equal to its number of holes.
    
    # Let h1 and h2 be the number of holes in patch 1 and patch 2, respectively.
    # We are given: h1 + h2 = total_holes.
    
    # The total number of logical qubits, K, is the sum from each independent patch.
    # K = k_patch1 + k_patch2
    # K = h1 + h2
    # Therefore, the total number of logical qubits is simply the total number of holes.
    
    total_logical_qubits = total_holes
    
    print("The number of logical qubits in a single surface code patch is equal to the number of holes in it.")
    print(f"We have {num_patches} patches and a total of {total_holes} holes to distribute among them.")
    print("The total number of logical qubits is the sum of the qubits from each patch.")
    
    # We can consider any distribution of the two holes. For example, let's put one hole in each patch.
    holes_in_patch1 = 1
    holes_in_patch2 = 1
    
    qubits_in_patch1 = holes_in_patch1
    qubits_in_patch2 = holes_in_patch2
    
    print("\nTo show the final equation, let's consider the case with one hole per patch:")
    print(f"Qubits in Patch 1 = {qubits_in_patch1}")
    print(f"Qubits in Patch 2 = {qubits_in_patch2}")
    print("The final calculation is:")
    print(f"{qubits_in_patch1} + {qubits_in_patch2} = {total_logical_qubits}")
    
    # The problem asks for the maximum. Since the total is constant regardless of distribution
    # (e.g., 2 holes in one patch and 0 in the other still gives 2 + 0 = 2), the maximum is 2.
    print(f"\nThe maximum number of logical qubits that can be encoded is {total_logical_qubits}.")


solve_logical_qubits()