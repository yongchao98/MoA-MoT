def solve_geology_problem():
    """
    This function analyzes the geological transition from Rock A to Rock B and determines the plausible causes.

    1.  **Geological Interpretation:**
        - Rock A is a dark, fine-grained rock, suggesting deposition in a deep, low-energy marine environment (e.g., deep sea shale).
        - Rock B is a lighter, coarser-grained rock, suggesting deposition in a shallower, higher-energy marine environment (e.g., nearshore sandstone).
        - The sequence from A to B represents a change from deeper to shallower water, which is called a marine regression.
        - A regression is caused by a relative fall in sea level.

    2.  **Evaluating the Options:**
        - We need to identify which events cause a fall in sea level.

        - **i. Onset of glaciation:** Glaciers form by trapping ocean water on land as ice, causing a global (eustatic) sea-level fall. This causes a regression. (Plausible)

        - **ii. Coastal subsidence:** Subsidence is the sinking of land. This causes a relative sea-level rise (transgression). (Not plausible)

        - **iii. Global warming:** Melting ice caps and thermal expansion of water cause a global sea-level rise (transgression). (Not plausible)

        - **iv. Decreased Mid-Ocean Ridge Activity:** Slower seafloor spreading causes the mid-ocean ridges to cool and contract, increasing the volume of ocean basins and causing a global (eustatic) sea-level fall. This causes a regression. (Plausible)

    3.  **Conclusion:**
        - The plausible explanations are the events that cause a sea-level fall.
    """
    
    plausible_events = ['i', 'iv']
    
    # Print the final answer in the required format
    answer = ",".join(plausible_events)
    print(answer)

solve_geology_problem()