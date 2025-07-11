import sys

def analyze_invasive_species():
    """
    Analyzes a list of species to determine which has had the largest negative
    ecosystem impact as an invasive introduced into New Mexico.
    """
    species_info = {
        "A. Apis mellifera (European Honey Bee)": 
            "Introduced and widespread. It can outcompete native pollinators, creating a negative ecological impact. However, its vital role in pollinating agricultural crops makes its overall impact complex and not typically categorized as the most negatively impactful invasive species.",
        
        "B. Aedes aegypti (Yellow Fever Mosquito)":
            "An invasive species that has become established in parts of New Mexico. It is a primary vector for dangerous diseases such as dengue, chikungunya, Zika, and yellow fever. Its impact on public health is a significant negative ecosystem impact, making it a major threat.",

        "C. Lycorma delicatula (Spotted Lanternfly)":
            "A highly destructive pest of many plants. As of recent data, it is a major problem in the eastern U.S. but is not yet established in New Mexico. Therefore, its current impact in New Mexico is minimal to none.",
        
        "D. Bombus pascuorum (Common Carder Bee)":
            "A species of bumblebee native to Europe. It is not known to be an established or impactful invasive species in North America, including New Mexico.",
            
        "E. Leptinotarsa decemlineata (Colorado Potato Beetle)":
            "Despite its name, this beetle is native to North America, likely originating in the region of Mexico and the Southwestern United States (including New Mexico). As a native species, it is not considered an 'introduced' invasive here.",
            
        "F. Maruca vitrata (Bean Pod Borer)":
            "A pest of legume crops found in tropical and subtropical regions. While it can cause agricultural damage, it is not considered the most significant invasive species in New Mexico in terms of broad ecosystem or public health impact."
    }
    
    print("Analysis of Invasive Species Impact in New Mexico:")
    print("="*50)
    for species, description in species_info.items():
        print(f"[*] {species}:\n    {description}\n")
        
    print("Conclusion:")
    print("Comparing the options, Aedes aegypti has the most severe negative impact in New Mexico. Its role as a vector for multiple debilitating human diseases presents a direct and significant threat to the ecosystem, including its human population.")

# Run the analysis
analyze_invasive_species()