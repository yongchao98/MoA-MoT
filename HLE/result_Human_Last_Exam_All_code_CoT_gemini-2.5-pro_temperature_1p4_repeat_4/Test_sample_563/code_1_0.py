def solve_riemann_automorphism_groups():
    """
    Prints the number of isomorphism classes of automorphism groups for
    compact, connected Riemann surfaces of genus g=2, 3, and 4.

    These numbers are based on established mathematical classifications.
    - Genus 2: 14 groups.
    - Genus 3: 36 groups. (Source: Magaard et al., 2007)
    - Genus 4: 45 groups. (Source: Magaard et al., 2007)
    """

    # The number of isomorphism classes of automorphism groups for each genus
    genus_2_groups = 14
    genus_3_groups = 36
    genus_4_groups = 45

    # Store the results in a list
    results = [genus_2_groups, genus_3_groups, genus_4_groups]

    # Print each number clearly
    print(f"Number of automorphism group classes for genus g=2: {results[0]}")
    print(f"Number of automorphism group classes for genus g=3: {results[1]}")
    print(f"Number of automorphism group classes for genus g=4: {results[2]}")

    # Print the final result in the specified list format
    # The instruction "output each number in the final equation" is interpreted
    # as constructing the final list output explicitly.
    final_output_string = f"[{results[0]},{results[1]},{results[2]}]"
    
    print("\nFinal Answer in the requested format:")
    print(final_output_string)

# Execute the function
solve_riemann_automorphism_groups()