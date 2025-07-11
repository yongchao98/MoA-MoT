def solve_wellington_mcq():
    """
    This function analyzes the historical statements about the Duke of Wellington
    and prints the numbers of the correct statements.
    """
    
    # Analysis of each statement based on historical facts.
    # A statement is considered correct if it aligns with historical records.
    
    # Statement 1: Correct. Wellington's logistical experience in India was foundational for his success in the Peninsular War.
    # Statement 2: Incorrect. His intelligence systems were famously successful in the Peninsular War (Europe).
    # Statement 3: Incorrect. This is an oversimplification. The influential divisional system was perfected in the Peninsula, not a direct transfer from the Deccan.
    # Statement 4: Incorrect. Sandhurst was founded much earlier, not by Wellington in 1829.
    # Statement 5: Incorrect. A strong negative claim that is false; the logistical lessons were highly significant.
    # Statement 6: Correct. The integration of local auxiliary forces was a hallmark of his command and became standard British imperial practice.
    # Statement 7: Incorrect. The opposite is true; his Indian experience was crucial for his commissariat reforms.
    # Statement 8: Correct. The use of 'flying columns' is a clear tactical concept transferred from India to other theaters like the Peninsula and Burma.
    # Statement 9: Incorrect. The 1813 Charter Act was driven by political and economic factors, not directly by Wellington's military principles.
    # Statement 10: Incorrect. His Indian civil administration experience was vital for managing territories and populations in the Peninsular War.
    
    correct_options = [1, 6, 8]
    
    # Sort the numbers and prepare the output string
    correct_options.sort()
    
    # We need to print the final equation. Since there is no equation, we will print the final answer in the format they wanted.
    answer_string = ",".join(map(str, correct_options))
    
    print("The final answer is:")
    print(f"<<<{answer_string}>>>")

solve_wellington_mcq()