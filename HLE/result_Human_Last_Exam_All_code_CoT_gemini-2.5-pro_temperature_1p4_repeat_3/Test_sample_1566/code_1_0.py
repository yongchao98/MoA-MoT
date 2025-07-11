# This script calculates the number of topologically distinct continua
# that satisfy the two given properties based on mathematical reasoning.

def solve():
    """
    Solves the problem by analyzing the topological properties.
    """

    # Let N be the number of topologically distinct continua.

    # Step 1: Analyze Property (1) and the definition of an end point.
    # The definition of an "end point" implies that the continuum X is "arc-like".
    # An arc-like continuum has at most 2 end points.
    # Property (1) states X has n end points where 1 < n < infinity.
    # Therefore, the number of end points must be exactly 2.
    num_endpoints = 2

    # Step 2: Analyze Property (2).
    # X has exactly 2 orbits under its homeomorphism group.
    # A known theorem in topology states that a continuum with a finite
    # number of orbits (like 2) must be locally connected.
    
    # Step 3: Combine the analyses.
    # X must be an arc-like, locally connected continuum with 2 end points.
    # The only such continuum, up to homeomorphism, is the closed interval [0, 1].
    
    # Step 4: Verify the solution.
    # The interval [0, 1] satisfies the properties:
    # 1. It has 2 end points ({0, 1}).
    # 2. It has 2 orbits under its homeomorphism group: the set of end points {0, 1}
    #    and the set of interior points (0, 1).
    
    # Step 5: Conclude the number of topologically distinct types.
    # Since any such continuum is homeomorphic to the interval [0, 1], there is
    # only one such topological type.
    N = 1
    
    # The "equation" is simply the final count.
    print(f"The number of topologically distinct continua is given by the final count.")
    print(f"Final Count: {N}")

solve()