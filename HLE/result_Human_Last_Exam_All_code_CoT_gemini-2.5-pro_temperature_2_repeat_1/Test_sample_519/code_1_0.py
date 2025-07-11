import sys

def solve():
    """
    This function provides the analysis for the three given categories fibered in groupoids.
    
    For X1, the Hilbert scheme of 11 points in A^3:
    - It's a Scheme (S).
    - It is separated (s).
    - It's known to be reducible for this degree and dimension.
    - Its dimension is a known result from the literature, dim = 35.
    
    For X2, the quotient stack [ (A^4 \\ V(xy-zw)) / C* ]:
    - Stabilizers are trivial, so it's an algebraic space, which is a DM-stack (DM).
    - It is separated (s).
    - The space being quotiented is irreducible, so the stack is irreducible (irr).
    - Dimension is dim(space) - dim(group) = 4 - 1 = 3.
    
    For X3, the Picard stack of a genus 7 curve:
    - Stabilizers are C*, which are not finite, so it is an Algebraic stack (A).
    - It is separated (s).
    - It's an infinite disjoint union of components (by degree), so it's not irreducible and not universally closed.
    - The dimension of each component is the genus of the curve, dim = 7.
    """
    
    # Properties for X1
    X1_type = "S"
    X1_props = ["s"]
    X1_dim = 35
    
    # Properties for X2
    X2_type = "DM"
    X2_props = ["s", "irr"]
    X2_dim = 3
    
    # Properties for X3
    X3_type = "A"
    X3_props = ["s"]
    X3_dim = 7
    
    # Format the output strings
    x1_str = f"[{X1_type}, {', '.join(X1_props)}, {X1_dim}]"
    # Python's list-to-string conversion has spaces after commas, 
    # need to remove them for the desired format.
    x1_str = x1_str.replace(" ", "")

    x2_str = f"[{X2_type}, {', '.join(X2_props)}, {X2_dim}]"
    x2_str = x2_str.replace(" ", "")
    
    x3_str = f"[{X3_type}, {', '.join(X3_props)}, {X3_dim}]"
    x3_str = x3_str.replace(" ", "")

    # In the problem description e.g. `[DM,s,7]` commas have no space after them
    # while f-string with join adds them. So let's re-format properly
    profile1 = f"[{X1_type},{','.join(X1_props)},{X1_dim}]"
    profile2 = f"[{X2_type},{','.join(X2_props)},{X2_dim}]"
    profile3 = f"[{X3_type},{','.join(X3_props)},{X3_dim}]"
    
    final_answer = f"{profile1} {profile2} {profile3}"
    print(final_answer)

solve()