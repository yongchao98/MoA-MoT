import pandas as pd

def find_most_impactful_invasive():
    """
    Analyzes a list of species to determine which has the largest negative
    ecosystem impact as an invasive in New Mexico.
    """

    # Data based on ecological records for New Mexico.
    # Impact Score is a simplified metric for this problem:
    # 0: Not an established invasive / No significant impact in NM.
    # 1: Minor impact (e.g., competition with native species) or is a native pest.
    # 2: Moderate ecosystem/agricultural impact.
    # 3: Major impact (e.g., disease vector, major ecological disruption).
    species_data = {
        'A': {'species': 'Apis mellifera', 'common_name': 'Honey Bee', 'impact_score': 1, 'notes': 'Introduced, competes with native pollinators but also beneficial for agriculture.'},
        'B': {'species': 'Aedes aegypti', 'common_name': 'Yellow Fever Mosquito', 'impact_score': 3, 'notes': 'Established invasive disease vector (Zika, dengue, etc.). Significant health and ecosystem impact.'},
        'C': {'species': 'Lycorma delicatula', 'common_name': 'Spotted Lanternfly', 'impact_score': 0, 'notes': 'Not established in New Mexico. A major threat, but not a current source of impact there.'},
        'D': {'species': 'Bombus pascuorum', 'common_name': 'Common Carder Bee', 'impact_score': 0, 'notes': 'Not a known significant invasive species in New Mexico.'},
        'E': {'species': 'Leptinotarsa decemlineata', 'common_name': 'Colorado Potato Beetle', 'impact_score': 1, 'notes': 'Native to the region (North America), not an "introduced" invasive. Is a major agricultural pest.'},
        'F': {'species': 'Maruca vitrata', 'common_name': 'Spotted Pod Borer', 'impact_score': 0, 'notes': 'Not considered a major invasive threat in New Mexico.'},
    }

    # Convert to a DataFrame for easy viewing
    df = pd.DataFrame.from_dict(species_data, orient='index')

    print("Analyzing species impact in New Mexico:")
    print("--------------------------------------------------")
    for choice, data in df.iterrows():
        print(f"Choice {choice}: {data['species']} -> Impact Score: {data['impact_score']}")
        print(f"   Note: {data['notes']}\n")


    # Find the species with the maximum impact score
    most_impactful_choice = df['impact_score'].idxmax()
    most_impactful_species = df.loc[most_impactful_choice]

    print("--------------------------------------------------")
    print("Conclusion:")
    print(f"The species with the highest impact score ({most_impactful_species['impact_score']}) in New Mexico is:")
    print(f"'{most_impactful_species['species']}' ({most_impactful_species['common_name']}).")
    print(f"Therefore, the correct choice is {most_impactful_choice}.")

if __name__ == '__main__':
    find_most_impactful_invasive()