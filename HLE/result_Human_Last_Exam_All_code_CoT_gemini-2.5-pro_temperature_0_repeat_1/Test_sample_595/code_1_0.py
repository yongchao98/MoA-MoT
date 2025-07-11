import math

def solve():
    """
    Finds the maximum number of grid squares the triangle's perimeter can cross.
    """
    max_k = 0
    best_u = 0

    # We only need to check for integer values of u from 1 to 17, as the floor
    # function is constant between integers. The case u=0 or u=18 corresponds
    # to aligning the triangle with the axes, which gives k=72.
    for u_int in range(1, 18):
        u = float(u_int)
        v_squared = 18**2 - u**2
        if v_squared < 0:
            continue
        v = math.sqrt(v_squared)

        # The vertices of the integer-grid triangle, f(P), where P are the
        # triangle's actual vertices. We assume one vertex is near the origin.
        pA_f = (0, 0)
        pB_f = (math.floor(u), math.floor(v))
        pC_f = (math.floor(-v), math.floor(u))

        # Calculate the L1 perimeter of the f(P) triangle
        k_AB = abs(pB_f[0] - pA_f[0]) + abs(pB_f[1] - pA_f[1])
        k_BC = abs(pC_f[0] - pB_f[0]) + abs(pC_f[1] - pB_f[1])
        k_CA = abs(pA_f[0] - pC_f[0]) + abs(pA_f[1] - pC_f[1])
        
        current_k = k_AB + k_BC + k_CA

        if current_k > max_k:
            max_k = current_k
            best_u = u

    # Now, print the calculation for the best case found
    u = best_u
    v = math.sqrt(18**2 - u**2)
    
    pA_f = (0, 0)
    pB_f = (math.floor(u), math.floor(v))
    pC_f = (math.floor(-v), math.floor(u))

    k_AB = abs(pB_f[0] - pA_f[0]) + abs(pB_f[1] - pA_f[1])
    k_BC = abs(pC_f[0] - pB_f[0]) + abs(pC_f[1] - pB_f[1])
    k_CA = abs(pA_f[0] - pC_f[0]) + abs(pA_f[1] - pC_f[1])

    print(f"The maximum number of squares crossed is k = {max_k}.")
    print("This can be achieved by orienting the triangle's legs with vectors (u, v) and (-v, u).")
    print(f"A maximal value is found for u = {int(u)}.")
    print(f"Calculation for u = {int(u)}:")
    print(f"v = sqrt(18^2 - {int(u)}^2) = sqrt({18**2 - int(u)**2}) ≈ {v:.2f}")
    print("\nThe integer-grid vertices are:")
    print(f"f(A) = {pA_f}")
    print(f"f(B) = (floor({u}), floor({v:.2f})) = {pB_f}")
    print(f"f(C) = (floor({-v:.2f}), floor({u})) = {pC_f}")
    print("\nNumber of squares crossed by each side:")
    print(f"Side AB: |{pB_f[0]} - {pA_f[0]}| + |{pB_f[1]} - {pA_f[1]}| = {abs(pB_f[0] - pA_f[0])} + {abs(pB_f[1] - pA_f[1])} = {k_AB}")
    print(f"Side BC: |{pC_f[0]} - {pB_f[0]}| + |{pC_f[1]} - {pB_f[1]}| = {abs(pC_f[0] - pB_f[0])} + {abs(pC_f[1] - pB_f[1])} = {k_BC}")
    print(f"Side CA: |{pA_f[0]} - {pC_f[0]}| + |{pA_f[1]} - {pC_f[1]}| = {abs(pA_f[0] - pC_f[0])} + {abs(pA_f[1] - pC_f[1])} = {k_CA}")
    print(f"\nTotal squares crossed k = {k_AB} + {k_BC} + {k_CA} = {max_k}")

solve()