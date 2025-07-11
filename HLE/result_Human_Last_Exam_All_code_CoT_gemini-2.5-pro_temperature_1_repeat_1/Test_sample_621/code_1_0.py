import sys

def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a rotation of A.
    B = A[i:] + A[:i]

    The algorithm has O(n) time complexity, which is optimal.
    """
    if len(A) != len(B) or not A:
        print("Input lists must be non-empty and of the same size.")
        return

    print(f"Given list A: {A}")
    print(f"Given list B: {B}")
    print("-" * 30)
    
    # --- Algorithm Explanation ---
    # The key insight is that the first element of the rotated list B,
    # which is B[0], must be the element at index 'i' in the original list A.
    # So, we have the relationship: B[0] == A[i].
    # To find 'i', we simply need to find the index of the element B[0] in list A.
    
    element_to_find = B[0]
    
    print("Step 1: Get the first element of list B.")
    print(f"   B[0] = {element_to_find}")
    
    print("\nStep 2: Find the index 'i' of this element in list A.")
    print(f"   The equation to solve is: A[i] = {element_to_find}")

    try:
        # The .index() method performs a linear search, which is O(n).
        found_i = A.index(element_to_find)
        
        print(f"\nSearching for {element_to_find} in A, we find it at index {found_i}.")
        print("-" * 30)
        print(f"Result: The rotation index is i = {found_i}")

        # --- Verification ---
        print("\nVerification:")
        reconstructed_B = A[found_i:] + A[:found_i]
        print(f"  Rotating A by i={found_i} gives us: {reconstructed_B}")
        print(f"  The original B was:              {B}")
        if reconstructed_B == B:
            print("  Success: The reconstructed B matches the original B.")
        else:
            print("  Failure: The lists do not match.")

    except ValueError:
        print(f"\nError: Element {element_to_find} (from B) was not found in A.")
        print("This implies B is not a valid rotation of A.")

if __name__ == '__main__':
    # Example Usage
    A = [3, 8, 1, 9, 4, 2, 7, 5]
    # Let's choose a rotation index to create B
    i_actual = 4
    B = A[i_actual:] + A[:i_actual]
    
    find_rotation_index(A, B)
