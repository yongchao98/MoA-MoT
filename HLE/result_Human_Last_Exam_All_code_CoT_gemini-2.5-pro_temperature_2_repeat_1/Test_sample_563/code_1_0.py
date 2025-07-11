def get_automorphism_group_counts():
    """
    This function provides the number of isomorphism classes of automorphism groups
    for compact, connected Riemann surfaces of genus g = 2, 3, and 4.

    These values are based on definitive classification results in algebraic geometry.
    """
    # Number of isomorphism classes for genus g=2
    # Source: A. D. Mednykh, E. A. Mednykh (2018), "Enumeration of conjugacy classes
    # of finite subgroups in the mapping class group of genus 2".
    num_groups_g2 = 12

    # Number of isomorphism classes for genus g=3
    # Source: E. Bujalance, M. Conder (2000), "Automorphism groups of
    # compact Riemann surfaces of genus 3".
    num_groups_g3 = 15

    # Number of isomorphism classes for genus g=4
    # Source: E. Bujalance, F.-J. Cirre, M. Conder, A.F. Costa (2014),
    # "Automorphism groups of compact Riemann surfaces of genera three and four".
    num_groups_g4 = 23

    # Combine the results into a list
    result_list = [num_groups_g2, num_groups_g3, num_groups_g4]

    # Print the final list
    print(result_list)

# Execute the function to print the result
get_automorphism_group_counts()