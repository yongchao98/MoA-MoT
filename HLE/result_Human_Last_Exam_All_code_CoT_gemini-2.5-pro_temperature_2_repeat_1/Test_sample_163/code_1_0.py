def determine_surveillance_protocol():
    """
    This function outlines the reasoning for selecting the best post-procedure
    surveillance program for a patient after SFA stenting.
    """

    # Step 1: Define the primary goal of surveillance.
    # The main goal after SFA stenting is to monitor for in-stent restenosis (the artery re-narrowing inside the stent).
    # This process, driven by neointimal hyperplasia, is most common within the first year.

    # Step 2: Compare surveillance modalities (ABI vs. Duplex Ultrasound).
    # - Ankle-Brachial Index (ABI): Measures overall blood pressure differences between the arm and ankle. While useful for initial diagnosis and tracking major hemodynamic changes, it is relatively insensitive for detecting early or moderate in-stent restenosis. A significant blockage can form before the ABI drops.
    # - Arterial Duplex Ultrasound: Combines imaging with Doppler flow measurement. It is the gold standard for stent surveillance because it allows for direct visualization of the stent and can precisely measure blood flow velocities to detect and quantify stenosis before it becomes severe or symptomatic.
    # - Conclusion: Surveillance programs based on duplex ultrasound are clinically superior. This points away from options A and C.

    # Step 3: Evaluate surveillance timing.
    # - Restenosis risk is highest in the first year. Therefore, surveillance should be more frequent during this period.
    # - Option E (annual visits) is too infrequent for early detection in the high-risk first year.
    # - Comparing options B and D:
    #   - B: Duplex at 1, 3, and 12 months. This schedule leaves a long 9-month gap between follow-ups during a critical period.
    #   - D: Duplex at 3, 6, 12 months, and 2 years. This provides more regular monitoring throughout the first year, including the 6-month mark which is a key time point for restenosis development. This schedule aligns well with consensus guidelines.

    # Step 4: Final Selection.
    # Option D combines the best imaging tool (arterial duplex) with the most clinically appropriate follow-up schedule
    # to effectively monitor for the primary complication of SFA stenting.

    print("Explanation of the chosen answer:")
    print("The most appropriate surveillance program uses the best tool at the right frequency.")
    print("1. Tool: Arterial duplex ultrasound is superior to ABI for detecting in-stent restenosis accurately and early.")
    print("2. Frequency: Follow-up should be most frequent in the first year. A schedule of 3, 6, and 12 months effectively monitors for the development of restenosis.")
    print("Choice D combines the correct tool and a robust, standard-of-care schedule.")


determine_surveillance_protocol()