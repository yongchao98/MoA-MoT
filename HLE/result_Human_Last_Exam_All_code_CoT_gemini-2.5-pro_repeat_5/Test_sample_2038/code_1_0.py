# First, you need to install the necessary libraries. You can do this by running:
# pip install spherogram
# spherogram also depends on snappy, which can be installed with:
# pip install snappy

import spherogram
import sys

def solve_knot_problem():
    """
    This script finds the number of knots with crossing number up to 13
    that admit two disjoint non-parallel embedded minimal genus Seifert surfaces.
    
    The condition for this is |det(K)| = 2*g(K) + 1.
    
    Note: This script may take a significant amount of time to run (potentially 30 minutes or more)
    due to the large number of knots with 12 and 13 crossings.
    """
    print("Starting the search for knots satisfying the condition |det(K)| = 2*g(K) + 1...")
    
    # This dictionary stores the number of prime knots for each crossing number.
    # These are standard values from knot theory databases (e.g., KnotInfo).
    knot_counts = {
        3: 1, 4: 1, 5: 2, 6: 3, 7: 7, 8: 21, 9: 49,
        10: 165, 11: 552, 12: 2176, 13: 9988
    }

    total_found = 0
    max_crossing = 13

    # Iterate through each crossing number from 3 to 13.
    for c in range(3, max_crossing + 1):
        num_knots_for_c = knot_counts.get(c, 0)
        if num_knots_for_c == 0:
            continue
            
        print(f"\n--- Checking {num_knots_for_c} knots with {c} crossings ---")
        found_at_c = 0
        
        # Iterate through each knot for the given crossing number.
        for i in range(1, num_knots_for_c + 1):
            knot_name = f'{c}_{i}'
            try:
                # Create the knot object.
                K = spherogram.Knot(knot_name)
                
                # Calculate the knot's genus.
                genus = K.genus()
                
                # Calculate the knot's determinant.
                determinant = K.determinant()
                
                # Check if the knot satisfies the Hirasawa-Sakuma condition.
                if determinant == 2 * genus + 1:
                    total_found += 1
                    found_at_c += 1
                    # As requested, output the numbers in the final equation.
                    print(f"Match found! Knot {knot_name}: det(K) = {determinant}, g(K) = {genus}. "
                          f"Equation: {determinant} = 2 * {genus} + 1.")

            except Exception as e:
                print(f"Could not process knot {knot_name}: {e}", file=sys.stderr)
        
        print(f"Found {found_at_c} matching knots with {c} crossings.")

    print("\n----------------------------------------------------")
    print("Search complete.")
    print(f"Total number of such knots with crossing number at most 13 is: {total_found}")
    
    # Returning the final answer for programmatic access.
    return total_found

if __name__ == '__main__':
    solve_knot_problem()