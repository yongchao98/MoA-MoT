import sys

def simulate_experiment():
    """
    This script simulates the expected results of a flow cytometry experiment
    to illustrate the importance of the anti-FLAG control as described in the problem.
    The data demonstrates the principle outlined in the correct answer.
    """

    # --- Simulated Experimental Data (Mean Fluorescence Intensity - MFI) ---

    # In the PBS (control) condition, both the anti-FLAG antibody (measuring total MUC1)
    # and the anti-glycan-MUC1 antibody (measuring the specific interaction)
    # should show a strong signal.
    pbs_mfi = {'anti_flag': 1215.5, 'anti_glycan_muc1': 1150.3}

    # In the 500 mM GalNAc (inhibitor) condition, we expect to see that the
    # total MUC1 level (measured by anti-FLAG) remains high and unchanged,
    # while the binding of the specific anti-glycan-MUC1 antibody is inhibited.
    galnac_mfi = {'anti_flag': 1205.2, 'anti_glycan_muc1': 180.8}

    # --- Output the results and interpretation ---

    print("Simulated Flow Cytometry Results (Mean Fluorescence Intensity)")
    print("="*75)
    print(f"{'Condition':<20} | {'Anti-FLAG MFI':<25} | {'Anti-Glycan-MUC1 MFI':<25}")
    print("-" * 75)
    
    # Print PBS control results
    print(f"{'PBS Control':<20} | {pbs_mfi['anti_flag']:<25.1f} | {pbs_mfi['anti_glycan_muc1']:<25.1f}")
    
    # Print GalNAc inhibitor results
    print(f"{'500 mM GalNAc':<20} | {galnac_mfi['anti_flag']:<25.1f} | {galnac_mfi['anti_glycan_muc1']:<25.1f}")
    
    print("="*75)
    print("\nInterpretation:")
    print("1. The 'Anti-FLAG MFI' is high and nearly identical in both the PBS and GalNAc conditions.")
    print("   This is the critical control. It confirms that the high concentration of GalNAc has NOT")
    print("   caused a decrease in the amount of MUC1 protein expressed on the cell surface.")
    
    print("\n2. The 'Anti-Glycan-MUC1 MFI' drops dramatically in the presence of GalNAc.")
    print("   Because the anti-FLAG control shows the MUC1 protein is still present, we can confidently")
    print("   conclude that this drop in signal is due to competitive inhibition of the antibody.")

    print("\nThis simulation demonstrates why it is essential to include the anti-FLAG antibody with the primary antibodies: to verify GalNAc has not altered surface expression of MUC1.")

if __name__ == '__main__':
    simulate_experiment()