def print_controller_formula():
    """
    Prints the derived formula for the set of all stabilizing controllers H_2(s).
    """
    
    # The Youla-Kucera parameter, K(s), can be any stable, proper rational function.
    # The derived controller H_2(s) stabilizes the plant H_1(s) = s/(s^2 - 1).
    
    nominator = "(s^2 - 1)K(s) + 5s + 4"
    denominator = "-s*K(s) - 4"
    
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("")
    print(f"               {nominator}")
    print(f"H_2(s)  =  -------------------------")
    print(f"               {denominator}")
    print("")
    print("where K(s) is any stable and proper rational function.")
    print("\nFor the controller H_2(s) to be proper, K(s) must also satisfy the condition:")
    print("lim_{s->inf} s*K(s) = -5")


print_controller_formula()

# The final answer in transfer function form is requested.
final_answer_nominator = "(s^2 - 1)K(s) + 5s + 4"
final_answer_denominator = "-(sK(s) + 4)" # Slightly rearranged for clarity
final_answer_string = f"({final_answer_nominator}) / ({final_answer_denominator})"

# Return the string representation for the answer format.
#<<<H_2(s) = ((s^2 - 1)K(s) + 5s + 4) / (-s*K(s) - 4)>>>