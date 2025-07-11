import sys

def solve_irritator_depressions():
    """
    This function calculates and prints the number of depressions
    in the lacrimal recess of Irritator challengeri based on
    published paleontological findings.
    """
    # According to the 2002 paper by Sues et al. on the skull of Irritator challengeri,
    # the lacrimal recess is described as having three main depressions.
    
    # We can represent each distinct depression.
    anterodorsal_depression = 1
    posterodorsal_depression = 1
    ventral_depression = 1
    
    # Calculate the total number of depressions.
    total_depressions = anterodorsal_depression + posterodorsal_depression + ventral_depression
    
    # Print the final result as a full sentence, showing each number in the equation.
    print(f"The lacrimal recess in Irritator challengeri contains {anterodorsal_depression} anterodorsal depression + {posterodorsal_depression} posterodorsal depression + {ventral_depression} ventral depression = {total_depressions} depressions in total.")

solve_irritator_depressions()