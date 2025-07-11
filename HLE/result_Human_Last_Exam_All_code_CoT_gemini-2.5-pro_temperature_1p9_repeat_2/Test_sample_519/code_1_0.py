def solve():
    """
    Calculates the properties of the three given categories fibered in groupoids.
    """

    # Properties for X1: Hilb^11(A^3)
    dim1 = 11 * 3
    props1 = ['S', 's', str(dim1)]
    
    # Properties for X2: [ (A^4 \ V(xy-zw)) / C* ]
    # dim = dim(space) - dim(group)
    dim2 = 4 - 1
    props2 = ['S', 's', 'irr', str(dim2)]

    # Properties for X3: Pic(C_0) for a genus 7 curve
    # dim = dim(H^1(C_0, O_C_0)) - dim(Aut(L)) = g - 1
    genus = 7
    dim3 = genus - 1
    props3 = ['A', 's', str(dim3)]

    # Format the output strings
    str1 = f"[{','.join(props1)}]"
    str2 = f"[{','.join(props2)}]"
    str3 = f"[{','.join(props3)}]"
    
    # Print the final result
    print(f"{str1} {str2} {str3}")

solve()
<<<[S,s,33] [S,s,irr,3] [A,s,6]>>>