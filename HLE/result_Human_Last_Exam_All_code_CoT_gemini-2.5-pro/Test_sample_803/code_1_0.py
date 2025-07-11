def main():
    """
    This script investigates the 'filled' property for a nonabelian group of order 54.
    The group is G = Q x_phi C_2, where Q is the nonabelian group UT(3,3) of order 27.
    
    According to theory, G is filled if and only if the set K = {u in Q | phi(u) = u^-1}
    is not a subgroup of Q. This is equivalent to K being a non-abelian set.
    This script verifies this condition for a specific choice of phi.
    """
    q = 3  # We are working with the field F_3 = {0, 1, 2}

    # --- 1. Define the group Q = UT(3,3) ---
    # Elements are tuples (x, y, z) representing upper unitriangular matrices.
    # Multiplication: (x1,y1,z1)*(x2,y2,z2) = (x1+x2, y1+y2, z1+z2+x1*y2) mod q
    
    Q_elements = []
    for x in range(q):
        for y in range(q):
            for z in range(q):
                Q_elements.append((x, y, z))

    # --- 2. Define group operations for Q ---
    def multiply(u1, u2):
        x1, y1, z1 = u1
        x2, y2, z2 = u2
        x_res = (x1 + x2) % q
        y_res = (y1 + y2) % q
        z_res = (z1 + z2 + x1 * y2) % q
        return (x_res, y_res, z_res)

    def inverse(u):
        x, y, z = u
        x_inv = (-x) % q
        y_inv = (-y) % q
        z_inv = (-z + x * y) % q
        return (x_inv, y_inv, z_inv)

    # --- 3. Define an automorphism phi of order 2 ---
    # We choose phi(x,y,z) = (-x, -y, z).
    def phi(u):
        x, y, z = u
        x_phi = (-x) % q
        y_phi = (-y) % q
        z_phi = z
        return (x_phi, y_phi, z_phi)

    # --- 4. Construct the set K = {u in Q | phi(u) = u^-1} ---
    K_set = []
    for u in Q_elements:
        if phi(u) == inverse(u):
            K_set.append(u)

    print(f"The group Q = UT(3, {q}) has {len(Q_elements)} elements.")
    print(f"The set K = {{u in Q | phi(u) = u^-1}} has {len(K_set)} elements.")
    print("The elements of K are:")
    # We sort for consistent output
    K_set.sort()
    print(K_set)

    # --- 5. Check if K is a subgroup by verifying closure ---
    # K is a subgroup if and only if it's closed under multiplication.
    is_subgroup = True
    counterexample = None
    for u1 in K_set:
        for u2 in K_set:
            prod = multiply(u1, u2)
            if prod not in K_set:
                is_subgroup = False
                counterexample = (u1, u2, prod)
                break
        if not is_subgroup:
            break

    # --- 6. Print the conclusion based on the subgroup property ---
    print("\n--- Checking if K is a subgroup ---")
    if is_subgroup:
        print("Result: K is closed under multiplication, so it is a subgroup.")
        print("Conclusion: The group G = Q x_phi C_2 is NOT a filled group.")
    else:
        u1, u2, prod = counterexample
        print("Result: K is not closed under multiplication, so it is not a subgroup.")
        print(f"Counterexample: u1 = {u1}, u2 = {u2}")
        print(f"u1 * u2 = {prod}, which is not in K.")
        phi_prod = phi(prod)
        inv_prod = inverse(prod)
        print(f"Verification: phi({prod}) = {phi_prod} and inverse({prod}) = {inv_prod}. They are not equal.")
        print("\nConclusion: Because K is not a subgroup, the group G = UT(3,3) x_phi C_2 is a FILLED group.")
        
    # --- 7. Alternative check: is K abelian? ---
    is_abelian = True
    non_commuting_pair = None
    for u1 in K_set:
        for u2 in K_set:
            prod12 = multiply(u1, u2)
            prod21 = multiply(u2, u1)
            if prod12 != prod21:
                is_abelian = False
                non_commuting_pair = (u1, u2)
                break
        if not is_abelian:
            break
            
    print("\n--- Alternative check: Is K abelian? ---")
    if is_abelian:
        print("Result: All elements in K commute. K is abelian.")
    else:
        u1, u2 = non_commuting_pair
        print("Result: K contains non-commuting elements. K is non-abelian.")
        print(f"Non-commuting pair in K: u1 = {u1}, u2 = {u2}")
        print(f"u1 * u2 = {multiply(u1, u2)}")
        print(f"u2 * u1 = {multiply(u2, u1)}")

if __name__ == '__main__':
    main()
