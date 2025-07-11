def check_homotopy_section_for_compact_manifold(manifold_name, euler_char, k):
    """
    For a compact manifold M, a homotopy section for the projection
    pi_{k, k+1}: conf_{k+1}(M) -> conf_k(M) exists if and only if
    the Euler characteristic chi(M) is equal to k.

    This function checks this condition for a given manifold and k.
    """
    
    print(f"Investigating the existence of a homotopy section for the map pi_{{{k}, {k+1}}} on the manifold M = {manifold_name}.")
    print(f"The manifold M is compact with Euler characteristic chi(M) = {euler_char}.")
    
    # The necessary and sufficient condition is chi(M) - k = 0.
    # We will check if this equation holds.
    result = euler_char - k
    
    print("\nThe condition to check is the equation: chi(M) - k = 0")
    print("Substituting the given values into the equation:")
    # The following print statement outputs each number in the final equation as requested.
    print(f"{euler_char} - {k} = {result}")
    
    print("\nConclusion:")
    if result == 0:
        print("The condition is satisfied. The map is expected to admit a homotopy section.")
    else:
        print("The condition is not satisfied. The map is not expected to admit a homotopy section.")


# Let's test the condition for the 2-sphere, M = S^2.
# The Euler characteristic of S^2 is 2.
sphere_euler_characteristic = 2

# We consider the map pi_{2,3}: conf_3(S^2) -> conf_2(S^2). Here, k = 2.
k_points = 2

check_homotopy_section_for_compact_manifold("S^2", sphere_euler_characteristic, k_points)
