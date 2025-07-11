import math

def solve_critical_correlation():
    """
    This function determines the 'critical amount of correlation' between
    input populations v and s required to balance potentiation and depression,
    allowing for selectivity to emerge.

    The balance condition is <dW/dt> = 0, which leads to the equality
    <r_i * v_k> = <r_i * s_k>.

    By expanding r_i in terms of its inputs and their statistical properties
    (mean activation 'mu', correlation 'C'), we derive the relationship:
    (mu - C) * (W_v - W_s) = 0.

    For selectivity to be possible, the network must support solutions where
    the weights to a neuron from the two populations are different (W_v != W_s).
    This implies that the other term in the equation must be zero, which gives
    the critical condition for the correlation.
    """

    # The inputs v and s are driven by a Poisson process with a given inter-event interval.
    # The average rate of activation 'mu' is the reciprocal of this interval.
    inter_event_interval = 150.0  # seconds

    # Calculate the average rate 'mu'
    mu = 1.0 / inter_event_interval

    # The critical condition for selectivity is mu - C = 0
    # Therefore, the critical amount of correlation C is equal to mu.
    C = mu

    # Print the explanation and the final result as per the instructions.
    print("The critical amount of correlation, C, is defined as the expectation value C = <v_k * s_k>.")
    print("The average rate of activation is denoted by \u03BC (mu), where \u03BC = <v_k> = <s_k>.")
    print("\nFrom the analysis, the condition that allows for weight selectivity (W^v \u2260 W^s) is found to be:")
    print("C = \u03BC")
    print("\nGiven the inter-event interval of {}s, the average rate \u03BC is 1 / {} = {:.6f} events/sec.".format(inter_event_interval, inter_event_interval, mu))
    
    print("\nTherefore, we can determine the numerical value for the critical correlation.")
    print("The final equation is C = \u03BC, where the values are:")
    print("C = {:.6f}".format(C))
    print("\u03BC = {:.6f}".format(mu))


solve_critical_correlation()

<<<0.006667>>>