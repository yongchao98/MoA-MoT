import pandas as pd

def calculate_simpsons_index(data):
    """
    Calculates Simpson's Diversity Index from a dictionary of species counts.
    Formula: D = 1 - (Σn(n-1)) / (N(N-1))
    
    Args:
        data (dict): A dictionary where keys are species names and values are their counts.
    """
    
    # Create a DataFrame for easy calculation and display
    df = pd.DataFrame(list(data.items()), columns=['Species', 'n (Count)'])
    
    # Calculate N, the total number of organisms
    N = df['n (Count)'].sum()
    
    # Calculate n-1 for each species
    df['n-1'] = df['n (Count)'] - 1
    
    # Calculate n(n-1) for each species
    df['n(n-1)'] = df['n (Count)'] * df['n-1']
    
    # Calculate Σn(n-1)
    sum_n_n_minus_1 = df['n(n-1)'].sum()
    
    print("Scenario: A sample is collected where only one species of bat is found.")
    print("Sample Data:")
    print(df[['Species', 'n (Count)']].to_string(index=False))
    print("-" * 30)
    
    print("Calculating the terms for the formula D = 1 - (Σn(n-1)) / (N(N-1))")
    
    # Avoid division by zero if N is 0 or 1
    if N <= 1:
        print("\nCannot calculate index because total individuals (N) is 1 or less.")
        # Simpson's index is technically undefined for N<2, but often considered 0 diversity.
        diversity_index = 0
    else:
        # Calculate N(N-1)
        N_N_minus_1 = N * (N - 1)
        
        print(f"Total individuals (N) = {N}")
        print(f"Σn(n-1) = {sum_n_n_minus_1}")
        print(f"N(N-1) = {N} * ({N}-1) = {N_N_minus_1}")
        print("-" * 30)

        # Calculate the final index
        diversity_index = 1 - (sum_n_n_minus_1 / N_N_minus_1)
        
        print("Final Calculation:")
        print(f"D = 1 - ({sum_n_n_minus_1} / {N_N_minus_1})")
        print(f"D = 1 - {sum_n_n_minus_1 / N_N_minus_1}")
        print(f"D = {diversity_index}")

    print("\nResult Analysis:")
    print("The mathematical calculation correctly yields an index of 0, confirming it's a valid mathematical outcome.")
    print("However, this result indicates a total lack of diversity (a monoculture) in the sample.")
    print("This contradicts the known ecological fact that several bat species exist on the island.")
    print("Therefore, the result is Mathematically valid, but not ecologically valid.")

# Hypothetical sample where the student only found 25 bats of a single species.
hypothetical_sample = {
    'Common Pipistrelle': 25
}

calculate_simpsons_index(hypothetical_sample)