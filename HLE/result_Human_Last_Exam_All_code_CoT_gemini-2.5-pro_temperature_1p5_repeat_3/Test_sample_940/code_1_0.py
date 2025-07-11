import sys

def evaluate_invasive_species():
    """
    Evaluates a list of species to determine which has had the largest
    negative impact as an introduced invasive in New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (Honey Bee)',
            'status': 'Introduced',
            'impact_in_nm': 'Widespread. While essential for agriculture, they are known to outcompete native pollinators for resources, which is a significant, albeit complex, negative ecosystem impact.'
        },
        'B': {
            'name': 'Aedes aegypti (Yellow Fever Mosquito)',
            'status': 'Introduced and Invasive',
            'impact_in_nm': 'Established in southern New Mexico. It is a vector for serious diseases like dengue, chikungunya, and Zika. Its presence poses a severe and direct threat to human and potentially animal health, representing a major negative impact.'
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted Lanternfly)',
            'status': 'Not established in NM',
            'impact_in_nm': 'A major threat, but it has not been found or established in New Mexico. Therefore, it has had no impact there yet.'
        },
        'D': {
            'name': 'Bombus pascuorum (Common Carder Bee)',
            'status': 'Not present in North America',
            'impact_in_nm': 'This species is native to Europe and is not found in New Mexico.'
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado Potato Beetle)',
            'status': 'Native',
            'impact_in_nm': 'This beetle is native to the region that includes New Mexico. While it is a major agricultural pest, it is not an "introduced" species in this context.'
        },
        'F': {
            'name': 'Maruca vitrata (Bean Pod Borer)',
            'status': 'Present, Agricultural Pest',
            'impact_in_nm': 'This moth is a pest on legume crops. Its impact is primarily agricultural and not considered as large of an ecosystem-wide negative impact as other options.'
        }
    }

    print("Analysis of Invasive Species Impact in New Mexico:")
    print("="*50)
    for choice, data in species_data.items():
        print(f"Option {choice}: {data['name']}")
        print(f"  - Status: {data['status']}")
        print(f"  - Impact: {data['impact_in_nm']}\n")

    print("Conclusion:")
    print("To be the correct answer, the species must be both introduced/invasive in New Mexico and have a large negative impact.")
    print("- Choices C, D, and E are eliminated because they are not established invasives in New Mexico (C, D) or are native (E).")
    print("- Between A (Apis mellifera) and B (Aedes aegypti), both are introduced and have negative impacts.")
    print("- However, the introduction of a potent disease vector like Aedes aegypti is generally considered a more severe and acutely negative impact than the resource competition from Apis mellifera.")
    print("\nTherefore, the species with the largest negative ecosystem impact from the list is B.")


if __name__ == '__main__':
    evaluate_invasive_species()