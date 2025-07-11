import textwrap

class AromaticCompound:
    """A simple class to represent the compounds and analyze their NMR properties."""
    def __init__(self, name, description, is_symmetric, thiophene_groups):
        """
        Initializes a compound.
        :param name: Name of the compound.
        :param description: A brief description.
        :param is_symmetric: Boolean, True if the molecule is symmetric.
        :param thiophene_groups: A list of strings describing the thiophene side groups.
                                 e.g., ['unbrominated', 'unbrominated']
        """
        self.name = name
        self.description = description
        self.is_symmetric = is_symmetric
        self.thiophene_groups = thiophene_groups

    def analyze_nmr(self):
        """Calculates the number of aromatic protons and expected NMR signals."""
        proton_definitions = {
            'unbrominated': ['C5-H (alpha)', 'C3-H (beta)'],
            'brominated_at_C5': ['C3-H (beta)']
        }

        total_protons = 0
        total_signals = 0
        
        # For asymmetric molecules, each proton on each group is unique
        if not self.is_symmetric:
            for group_type in self.thiophene_groups:
                protons_in_group = proton_definitions.get(group_type, [])
                total_protons += len(protons_in_group)
                total_signals += len(protons_in_group)
            return total_protons, total_signals

        # For symmetric molecules, count unique proton types
        # Assumes all groups are identical
        if self.is_symmetric and self.thiophene_groups:
            group_type = self.thiophene_groups[0]
            protons_per_group = len(proton_definitions.get(group_type, []))
            signals_per_group = protons_per_group
            num_groups = len(self.thiophene_groups)
            
            total_protons = protons_per_group * num_groups
            total_signals = signals_per_group

            return total_protons, total_signals
        
        return 0, 0

def solve_chemistry_problem():
    """Executes the reasoning to solve the user's problem."""
    
    # Define the molecules involved in the reaction
    sm_name = "Starting Material"
    sm_desc = ("2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]"
               "isoindole-4,6(5H)-dione")
    starting_material = AromaticCompound(sm_name, sm_desc, is_symmetric=True, 
                                         thiophene_groups=['unbrominated', 'unbrominated'])

    dibromo_product = AromaticCompound("Dibrominated Product", "Product of reaction with 2 eq. NBS",
                                       is_symmetric=True, 
                                       thiophene_groups=['brominated_at_C5', 'brominated_at_C5'])

    monobromo_product = AromaticCompound("Monobrominated Product", "Intermediate product",
                                         is_symmetric=False, 
                                         thiophene_groups=['unbrominated', 'brominated_at_C5'])
    
    # --- Print out the step-by-step thinking process ---
    
    print("Here is a step-by-step analysis to identify the unknown compound:\n")
    print("--- Step 1: Analyze the Starting Material and Reaction ---")
    print(f"The starting material is {starting_material.description}.")
    print("The core of this molecule has no protons. The reactive sites for electrophilic bromination "
          "with NBS are the C-H bonds on the two pendant thiophene rings.")
    print("On each thiophene ring, the most reactive position is the unsubstituted alpha-position (C5), "
          "followed by the beta-position (C3).\n")

    print("--- Step 2: Correlate Structure with NMR Data ---")
    print("The key to identifying the product is the ¹H-NMR data: 'three peaks that are larger than 6.0 ppm'.")
    print("Let's analyze the expected number of aromatic proton signals for each possible compound.\n")

    print("--- Step 3: Evaluate Each Possible Compound ---")

    # Analysis of Starting Material
    sm_protons, sm_signals = starting_material.analyze_nmr()
    print(f"A) {starting_material.name}:")
    print(f"   - Structure: It has two identical, unbrominated thiophene rings. The molecule is symmetric.")
    print(f"   - Aromatic Protons: It has a total of {sm_protons} aromatic protons (2 on each ring).")
    print(f"   - Expected NMR Signals: Due to symmetry, the protons on the two rings are equivalent, giving {sm_signals} distinct signals.")
    print(f"   - Verdict: Does not match the 'three peaks' observation.\n")

    # Analysis of Dibrominated Product
    di_protons, di_signals = dibromo_product.analyze_nmr()
    print(f"B) {dibromo_product.name}:")
    print(f"   - Structure: This is the expected product from reacting with 2 equivalents of NBS. Both C5 positions are brominated.")
    print(f"   - Aromatic Protons: The two C5 protons are replaced by bromine. Only the two C3 protons remain, for a total of {di_protons} aromatic protons.")
    print(f"   - Expected NMR Signals: Due to symmetry, the two remaining C3 protons are equivalent, giving only {di_signals} signal.")
    print(f"   - Verdict: Does not match the 'three peaks' observation.\n")

    # Analysis of Monobrominated Product
    mono_protons, mono_signals = monobromo_product.analyze_nmr()
    print(f"C) {monobromo_product.name}:")
    print(f"   - Structure: One thiophene ring is brominated at C5, while the other is unreacted. The molecule is now asymmetric.")
    print(f"   - Aromatic Protons (and their corresponding NMR signals):")
    print(f"     - Unreacted thiophene: 1 proton at C5' and 1 proton at C3'. These are unique. (2 signals)")
    print(f"     - Brominated thiophene: 1 proton at C3''. This is also unique. (1 signal)")
    print(f"   - Total Protons = 1 + 1 + 1 = {mono_protons} protons.")
    print(f"   - Expected NMR Signals = 2 + 1 = {mono_signals} distinct signals, each integrating to one proton.")
    print(f"   - Verdict: This perfectly matches the 'three peaks' observation.\n")
    
    print("--- Step 4: Conclusion ---")
    print("The observation that the reaction was sluggish with 2 eq. of NBS but proceeded after adding more suggests that the reaction kinetics are slow. By stopping the reaction after one hour, you likely isolated the monobrominated intermediate before it could fully convert to the dibrominated product.")
    print("The ¹H-NMR data definitively points to the asymmetrically substituted monobrominated compound.\n")
    
    print("=" * 70)
    print("FINAL ANSWER")
    print("=" * 70)
    print("The new spot isolated from the reaction is the monobrominated product.")
    
    product_name = ("2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(4-(2-ethylhexyl)thiophen-2-yl)"
                    "-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione")
    
    print("\nFull IUPAC Name:")
    print('\n'.join(textwrap.wrap(product_name, 70)))
    
    return product_name

# Run the analysis and get the final answer for the '<<<' tag.
final_answer = solve_chemistry_problem()

print(f"\n<<<{final_answer}>>>")
