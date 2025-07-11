def get_automorphism_group_counts():
    """
    This function provides the number of isomorphism classes of automorphism groups
    for compact, connected Riemann surfaces of genus g=2, 3, and 4.
    
    These numbers are based on established mathematical classification results, specifically
    the modern and comprehensive data compiled by Thomas Breuer.
    """
    
    # Number of isomorphism classes of automorphism groups for a Riemann surface of genus 2
    num_g2 = 13
    
    # Number of isomorphism classes of automorphism groups for a Riemann surface of genus 3
    num_g3 = 23
    
    # Number of isomorphism classes of automorphism groups for a Riemann surface of genus 4
    num_g4 = 38
    
    # The final list of counts for g=2, g=3, and g=4
    final_list = [num_g2, num_g3, num_g4]
    
    # Per instructions, outputting each number in the final equation.
    print(f"Number of groups for genus g=2 is {num_g2}")
    print(f"Number of groups for genus g=3 is {num_g3}")
    print(f"Number of groups for genus g=4 is {num_g4}")
    
    # Final answer in the specified list format
    print(final_list)

get_automorphism_group_counts()