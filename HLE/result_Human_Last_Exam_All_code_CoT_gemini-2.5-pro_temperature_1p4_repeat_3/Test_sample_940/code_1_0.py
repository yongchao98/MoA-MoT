import pandas as pd

def find_most_impactful_invasive():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive species in New Mexico.
    """
    data = {
        'Species': [
            'A. Apis mellifera',
            'B. Aedes aegypti',
            'C. Lycorma delicatula',
            'D. Bombus pascuorum',
            'E. Leptinotarsa decemlineata',
            'F. Maruca vitrata'
        ],
        'Common Name': [
            'European Honey Bee',
            'Yellow Fever Mosquito',
            'Spotted Lanternfly',
            'Common Carder Bee',
            'Colorado Potato Beetle',
            'Bean Pod Borer'
        ],
        'Is_Invasive_In_NM': [
            True,
            True,
            False,
            False,
            False,  # It is native to the region
            True
        ],
        'Impact_Score': [
            4,  # Competes with native pollinators, but also beneficial
            9,  # Major public health vector for serious diseases
            0,  # Not yet established in New Mexico
            0,  # Not established in North America
            0,  # Native to NM/CO, not an "introduced" invasive
            3   # Agricultural pest, but less impact than Aedes aegypti
        ],
        'Reasoning': [
            "Introduced and widespread. Competes with native pollinators but is also a critical agricultural pollinator. Its negative impact is debated.",
            "An invasive species in southern New Mexico that vectors serious diseases like Dengue, Zika, and Chikungunya, posing a major public health risk.",
            "A highly destructive pest in the eastern US, but it is not yet established in New Mexico.",
            "A European bumblebee that is not invasive in North America.",
            "This species is native to the region of Colorado and New Mexico. It is an agricultural pest but not considered an introduced invasive species here.",
            "An agricultural pest on legumes, but its overall negative impact in New Mexico is less severe than disease-vectoring species."
        ]
    }

    df = pd.DataFrame(data)

    # Filter for species that are actually invasive in New Mexico
    invasive_in_nm_df = df[df['Is_Invasive_In_NM'] == True]

    # Find the species with the highest impact score from the filtered list
    most_impactful = invasive_in_nm_df.loc[invasive_in_nm_df['Impact_Score'].idxmax()]

    print("Analysis of Invasive Species Impact in New Mexico:")
    print("-" * 50)
    for index, row in df.iterrows():
        print(f"Species: {row['Species']} ({row['Common Name']})")
        print(f"  - Status: {row['Reasoning']}\n")

    print("-" * 50)
    print("Conclusion:")
    print(f"The species with the largest negative impact as an invasive species in New Mexico is:")
    print(f"{most_impactful['Species']} ({most_impactful['Common Name']})")
    print(f"Reason: {most_impactful['Reasoning']}")

if __name__ == '__main__':
    find_most_impactful_invasive()
<<<B>>>