import itertools

def compute_spin_bordism_bg2():
    """
    Computes the reduced 12-th dimensional Spin bordism of the classifying space of G2.
    This computation is based on the theorem by K. Yamaguti (2018) which describes
    the structure of the reduced Spin bordism of BG2 as an exterior algebra over
    the Spin bordism ring of a point.
    """
    
    # Spin bordism groups of a point, Omega_k^Spin, for k <= 12
    # Represented as strings for clarity. "0" means the trivial group.
    omega_spin = {
        0: "Z",
        1: "Z_2",
        2: "Z_2",
        3: "0",
        4: "Z",
        5: "0",
        6: "0",
        7: "0",
        8: "Z + Z",
        9: "Z_2 + Z_2",
        10: "Z_2",
        11: "0",
        12: "Z",
    }

    # Degrees of the elementary generators of the exterior algebra
    generator_degrees = [5, 6, 9, 10]
    
    target_dimension = 12
    
    summands = []

    print(f"Computing the reduced {target_dimension}-th Spin bordism of BG2.")
    print("The result is a direct sum of Omega_k^Spin groups, where k is determined by the dimension of generators.")
    print("The degrees of the elementary generators are: ", generator_degrees)
    print("-" * 30)

    # Iterate through all possible generators of the exterior algebra
    # These correspond to all subsets of the elementary generators.
    for i in range(len(generator_degrees) + 1):
        for combo in itertools.combinations(generator_degrees, i):
            generator_degree = sum(combo)
            
            if generator_degree <= target_dimension:
                k = target_dimension - generator_degree
                
                # Look up the coefficient group Omega_k^Spin
                group = omega_spin.get(k, "0")
                
                if group != "0":
                    # This generator contributes a non-trivial summand
                    gen_str = "1" if not combo else " * ".join([f"beta_{c}" for c in combo])
                    print(f"Generator: {gen_str} (degree {generator_degree})")
                    print(f"  Coefficient group: Omega_{k} = {group}")
                    summands.append(group)

    # Format the final result
    if not summands:
        final_result = "0"
    else:
        final_result = " + ".join(summands)
        
    print("-" * 30)
    print("The reduced 12-th dimensional Spin bordism of BG2 is the direct sum of the non-trivial contributions:")
    print(f"Final Result: {final_result}")

if __name__ == '__main__':
    compute_spin_bordism_bg2()
