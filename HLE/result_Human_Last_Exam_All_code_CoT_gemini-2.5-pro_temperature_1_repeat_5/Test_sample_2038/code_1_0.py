import spherogram
import requests
from collections import defaultdict
import sys

def solve_knot_problem():
    """
    This script calculates the number of 2-bridge knots with crossing number at most 13
    that admit two disjoint non-parallel embedded minimal genus Seifert surfaces.
    """
    print("Starting the calculation. This may take a few minutes, especially the first time it's run...")

    # A knot has the desired property iff it's a fibered 2-bridge knot.
    # A knot is 2-bridge iff its determinant is odd.
    # All 2-bridge knots are alternating.

    # Part 1: Use spherogram for knots with crossing number <= 12
    counts_by_crossing = defaultdict(int)
    # The alternating_knots catalog is complete up to 12 crossings.
    try:
        knot_iterator = spherogram.alternating_knots(max_crossings=12)
        total_le12 = len(spherogram.alternating_knots.catalog_up_to(12))
        
        for i, K in enumerate(knot_iterator):
            # Provide progress update
            if i % 100 == 0:
                sys.stdout.write(f"\rProcessing knots up to 12 crossings: {i+1}/{total_le12}...")
                sys.stdout.flush()

            c = K.crossing_number()
            # Check if it's a 2-bridge knot (odd determinant) and if it's fibered.
            if K.determinant() % 2 != 0 and K.is_fibered():
                counts_by_crossing[c] += 1
        print(f"\rProcessing knots up to 12 crossings: Done.               ")

    except Exception as e:
        print(f"\nAn error occurred while using spherogram: {e}")
        print("Please ensure 'spherogram' and its dependencies are correctly installed.")
        return

    # Part 2: Fetch data from KnotInfo for 13-crossing knots
    print("Fetching data for 13-crossing knots from KnotInfo database...")
    count13 = 0
    url = 'http://knotinfo.math.indiana.edu/data/13_alternating.txt'
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.text.splitlines()

        header = data[0].split()
        det_idx = header.index('det')
        fib_idx = header.index('fibered')

        for line in data[1:]:
            fields = line.split()
            if not fields:
                continue
            
            determinant = int(fields[det_idx])
            fibered_status = fields[fib_idx]

            if determinant % 2 != 0 and fibered_status == 'F':
                count13 += 1
        
        counts_by_crossing[13] = count13
        print("Successfully processed data for 13-crossing knots.")

    except requests.exceptions.RequestException as e:
        print(f"\nError: Could not fetch data for 13-crossing knots.")
        print(f"Details: {e}")
        print("Cannot compute the final count. Please check your internet connection.")
        return
    except (ValueError, IndexError):
        print("\nError: Could not parse the data file from KnotInfo.")
        print("The data format may have changed. Cannot compute the final count.")
        return

    # Part 3: Print the results
    print("\n--- Results ---")
    equation_parts = []
    for c in sorted(counts_by_crossing.keys()):
        count = counts_by_crossing[c]
        print(f"Number of such knots with crossing number {c}: {count}")
        equation_parts.append(str(count))

    total_count = sum(counts_by_crossing.values())
    equation_str = " + ".join(equation_parts)
    
    print("\nFinal calculation:")
    print(f"{equation_str} = {total_count}")
    
    # Required final answer format
    print(f"\n<<<{total_count}>>>")


if __name__ == '__main__':
    solve_knot_problem()
