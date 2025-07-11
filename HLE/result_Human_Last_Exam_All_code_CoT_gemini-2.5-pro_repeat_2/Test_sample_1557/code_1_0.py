import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def solve_history_mcq():
    """
    This function evaluates the historical statements about the Duke of Wellington's career
    and identifies the correct ones based on established historical facts.
    """

    # Evaluation of each statement:
    # 1. Correct. Wellington's logistical experience in India, particularly with the mobile commissariat using bullock carts,
    #    was foundational to his success in the Peninsular War. Historians widely acknowledge this continuity.
    # 2. Incorrect. While he developed a superb intelligence network in India, he also created one of the most effective
    #    intelligence systems of the Napoleonic Wars in the Peninsula, disproving the claim it was "never successfully implemented in European theaters."
    # 3. Incorrect. The divisional system perfected in the Peninsula became the influential model, not a specific Deccan structure.
    #    The post-1815 army organization was influenced by the Peninsular experience, but the statement is an oversimplification.
    # 4. Incorrect. The Royal Military College, Sandhurst was founded in 1801, not 1829 by Wellington. The East India Company had its own seminary at Addiscombe for India-specific training.
    # 5. Incorrect. This is the opposite of the truth. Logistical lessons from India (e.g., use of local transport, long-distance supply)
    #    were highly relevant and applied in other colonial campaigns (e.g., Burma, Afghanistan).
    # 6. Correct. The practice of raising and integrating local troops (sepoys) under British command, which Wellington worked with extensively in India,
    #    became a template for British imperial forces across the globe, including in Africa and Southeast Asia.
    # 7. Incorrect. This is definitively false and contradicts historical consensus. His Indian service was a crucial learning ground for the logistical challenges he later solved in the Peninsula.
    # 8. Correct. The concept of fast-moving, independent "flying columns" was used by Wellington in India against the Marathas, adapted for the Peninsula,
    #    and was a key tactic in later colonial conflicts like the First Anglo-Burmese War (1824-1826).
    # 9. Incorrect. The East India Company Charter Act of 1813 was driven by economic and political pressures in Britain (e.g., free trade lobbyists),
    #    not by Wellington's military organizational principles.
    # 10. Incorrect. Wellington's methods of managing civilian populations—paying for supplies, preventing looting, and maintaining order—were developed in India
    #     and applied directly and successfully in the Peninsula, contributing significantly to his campaign's success.

    correct_options = [1, 6, 8]
    
    # Sort the options and format for output
    correct_options.sort()
    
    # The final print statement for the user.
    print(','.join(map(str, correct_options)))

solve_history_mcq()

# Restore the original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output = captured_output.getvalue()

# Final printing of the captured output
# This ensures the only printed output is the one from the function
print(output.strip())