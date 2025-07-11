import urllib.request
import re
import sys

def count_automorphism_groups_by_genus(genus):
    """
    Counts the number of isomorphism classes of automorphism groups for a given genus.

    This function fetches data from S. Allen Broughton's "Atlas of Riemann Surfaces Actions"
    (http://www.math.wku.edu/~broughtn/data/data.html), which lists the possible
    automorphism groups for Riemann surfaces of a given genus. Each group isomorphism
    class is uniquely identified by its GAP SmallGroup ID, e.g., '[24,12]' for S4.

    Args:
        genus (int): The genus of the Riemann surface (g >= 2).

    Returns:
        int: The number of unique isomorphism classes of automorphism groups.
             Returns None if data fetching or parsing fails.
    """
    url = f"http://www.math.wku.edu/~broughtn/data/genus{genus}.dat"
    try:
        with urllib.request.urlopen(url) as response:
            if response.status != 200:
                print(f"Error: Failed to fetch data for genus {genus}. HTTP Status: {response.status}", file=sys.stderr)
                return None
            data = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error: Could not retrieve or read data from {url}. Details: {e}", file=sys.stderr)
        return None

    # The data file format includes the GAP SmallGroup ID in the format '[order,index]'.
    # Example line: 4 3 24 2 2 3 : 24 10: [24,12] S4 : [1,4][1,1][2,3]
    # The regular expression finds the colon preceding the group ID, then captures
    # the ID itself. This is a reliable way to isolate the group identifiers.
    group_ids = re.findall(r":\s*\[(\d+,\d+)\]", data)

    # A set is used to store the group IDs, which automatically handles uniqueness.
    # The size of the set gives the number of isomorphism classes.
    unique_ids = set(group_ids)
    
    return len(unique_ids)

def main():
    """
    Main function to calculate the number of automorphism groups for
    genera 2, 3, and 4 and print the result.
    """
    genera = [2, 3, 4]
    results = []
    
    for g in genera:
        count = count_automorphism_groups_by_genus(g)
        if count is not None:
            results.append(count)
        else:
            # In case of an error, print an error message and exit.
            print(f"Could not compute the result for genus {g}.", file=sys.stderr)
            return

    # Print the final result in the specified format [num_g2, num_g3, num_g4]
    print(f"[{results[0]},{results[1]},{results[2]}]")

if __name__ == "__main__":
    main()
