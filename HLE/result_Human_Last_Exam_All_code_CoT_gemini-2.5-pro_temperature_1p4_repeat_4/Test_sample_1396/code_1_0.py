def calculate_upper_bound():
    """
    Calculates the upper bound on the number of cuts for a connected
    envy-free cake cutting allocation for 4 agents based on the
    Brânzei and Nisan (2022) algorithm.
    """
    
    # Number of agents in the problem
    num_agents = 4
    
    # Number of cuts each agent makes in the protocol
    cuts_per_agent = 2
    
    # Number of cuts the referee makes
    referee_cuts = 2
    
    # The upper bound O is the total number of cuts required in the worst case.
    # It's calculated as (number of agents * cuts per agent) + referee cuts.
    upper_bound_O = (num_agents * cuts_per_agent) + referee_cuts
    
    print("Based on the state-of-the-art algorithm by Brânzei and Nisan (2022) for 4 agents, the upper bound 'O' on the number of cuts is calculated.")
    print("The final equation for the total cuts is:")
    print(f"({num_agents} * {cuts_per_agent}) + {referee_cuts} = {upper_bound_O}")
    print(f"\nThe most realistic upper bound O is {upper_bound_O}.")

calculate_upper_bound()