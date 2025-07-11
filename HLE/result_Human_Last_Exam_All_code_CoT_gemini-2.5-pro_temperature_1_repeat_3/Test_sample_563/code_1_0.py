def get_automorphism_group_counts():
    """
    This function provides the number of isomorphism classes of automorphism groups
    for compact, connected Riemann surfaces of genus g = 2, 3, and 4.

    These numbers are established results from the mathematical literature on the
    classification of finite group actions on surfaces.
    """

    # For genus g=2, there are 12 isomorphism classes of automorphism groups.
    num_groups_genus_2 = 12

    # For genus g=3, there are 36 isomorphism classes of automorphism groups.
    # This was famously classified by A. Kuribayashi and K. Komiya.
    num_groups_genus_3 = 36

    # For genus g=4, there are 23 isomorphism classes of automorphism groups.
    # This result is based on classifications by S. A. Broughton and others.
    num_groups_genus_4 = 23

    # The final answer is a list containing these numbers.
    result = [num_groups_genus_2, num_groups_genus_3, num_groups_genus_4]

    # Print the final result in the required format.
    # The instruction "output each number in the final equation!" is interpreted
    # as printing the list containing each of the required numbers.
    print(result)

get_automorphism_group_counts()