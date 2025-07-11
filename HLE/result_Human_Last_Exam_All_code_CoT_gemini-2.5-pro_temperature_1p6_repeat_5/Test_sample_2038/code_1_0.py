import pandas as pd
import io
import requests

def solve_knot_problem():
    """
    This function calculates the number of 2-bridge knots with a crossing number
    at most 13 that have two disjoint non-parallel embedded minimal genus Seifert surfaces.
    
    This is equivalent to counting the number of non-alternating 2-bridge knots
    with a crossing number of at most 13.
    
    The function fetches data from the KnotInfo database, filters it according to the
    problem's criteria, and prints the detailed results.
    """
    
    # URL of the KnotInfo data CSV file
    url = "https://knotinfo.math.indiana.edu/knot_data.csv"

    print("Step 1: The knot property in question is equivalent to the knot being non-alternating.")
    print("Step 2: The plan is to count non-alternating 2-bridge knots with crossing number <= 13 by analyzing the KnotInfo database.\n")
    
    try:
        # Fetch the data from the URL
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an HTTPError for bad responses
        csv_data = response.content
        
        # Read the CSV data into a pandas DataFrame
        df = pd.read_csv(io.StringIO(csv_data.decode('utf-8')))

        # --- Data Cleaning and Filtering ---
        # Convert relevant columns to numeric types, coercing errors to NaN
        df['CrossingNum'] = pd.to_numeric(df['CrossingNum'], errors='coerce')
        df['BridgeNum'] = pd.to_numeric(df['BridgeNum'], errors='coerce')

        # Drop any rows where conversion failed for essential columns
        df.dropna(subset=['CrossingNum', 'BridgeNum', 'Alternating'], inplace=True)

        # Apply the filters to find the knots that match our criteria
        filtered_knots = df[
            (df['CrossingNum'] <= 13) &
            (df['BridgeNum'] == 2) &
            (df['Alternating'] == 'No')
        ]

        # --- Calculation and Output ---
        if filtered_knots.empty:
            print("No knots were found that match the criteria.")
            return
            
        # Group by crossing number and count the knots in each group
        counts_by_crossing = filtered_knots.groupby('CrossingNum').size()
        
        print("Step 3: Count of such knots found for each crossing number:")
        
        total_count = 0
        equation_parts = []
        
        # We start from crossing number 3, as it's the first non-trivial knot.
        for c in range(3, 14):
            # .get(c, 0) handles crossing numbers with no matching knots
            count_at_c = counts_by_crossing.get(float(c), 0)
            if count_at_c > 0:
                print(f"  - Crossing number {c}: {count_at_c}")
                total_count += count_at_c
                equation_parts.append(str(count_at_c))
        
        print("\nStep 4: The final count is the sum of the counts for each crossing number.")
        
        # Construct and print the final equation
        equation_str = " + ".join(equation_parts)
        print(f"Final calculation: {equation_str} = {total_count}")
        
    except requests.exceptions.RequestException as e:
        print(f"Error: Failed to fetch data from KnotInfo database. Please check your internet connection.\n{e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Execute the main function
solve_knot_problem()