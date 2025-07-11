import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def solve_epidemiological_model(params):
    """
    Solves the system of ODEs for the given parameters.
    """
    # Nominal Parameters
    p = {
        'a_i': 6, 'a_r': 21, 'f_s': 0.2, 'mu': 1/45625, 'mu_n': 1/2100,
        'mu_s': 1/600, 'mu_h': 1/2400, 'beta_h': 0.1, 'r_b': 0, 'c_h': 1,
        'c_l': 1, 'q_s': 60, 'q_l': 56, 'q_f': 0.5, # Assuming a nominal q_f
        'beta_n_val': 0.5, 'beta_s_val': 0.5
    }
    # Update with any varied parameters
    p.update(params)

    # ODE system function
    def model(t, y, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        # Total population
        T = max(0, S + E + In + Is + R)
        if T == 0:
            return np.zeros(10)

        # Time-dependent contact rates
        if p['q_s'] <= t <= p['q_s'] + p['q_l']:
            beta_n = p['q_f'] / 2
            beta_s = p['q_f'] / 2
        else:
            beta_n = p['beta_n_val']
            beta_s = p['beta_s_val']

        # Infection terms
        infection_from_Is = beta_s * S * max(0, Is - H) / T
        infection_from_H = p['beta_h'] * S * H / T
        infection_from_In = beta_n * S * In / T
        total_infection = infection_from_Is + infection_from_H + infection_from_In

        # Derivatives
        dSdt = -total_infection - p['mu'] * S
        dEdt = total_infection - (1/p['a_i'] + p['mu']) * E
        dIndt = E * (1 - p['f_s']) / p['a_i'] - In / p['a_r'] - p['mu_n'] * In
        
        dIsdt_infection_term = E * p['f_s'] / p['a_i']
        dIsdt = dIsdt_infection_term - Is / p['a_r'] - p['mu_h'] * H - p['mu_s'] * max(0, Is - H)
        
        # Hospitalization term logic
        if H < B:
            hospitalization = min(B - H, dIsdt_infection_term)
        else:
            hospitalization = 0

        dHdt = hospitalization - H / p['a_r'] - p['mu_h'] * H
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * max(0, Is - H) + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Initial conditions
    y0 = [999999, 0, 1, 0, 0, 0, 0, 2900000, 0, 0] # S, E, In, Is, H, R, D, B, Ch, Cl
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], 500)
    
    sol = solve_ivp(model, t_span, y0, args=(p,), t_eval=t_eval, method='RK45')
    
    # Recalculate T for output
    S, E, In, Is, H, R, D, B, Ch, Cl = sol.y
    T = np.maximum(0, S + E + In + Is + R)

    results = {
        'S': S, 'E': E, 'In': In, 'Is': Is, 'H': H, 'R': R, 'D': D, 
        'B': B, 'Ch': Ch, 'Cl': Cl, 'T': T, 't': sol.t
    }
    return results

def get_solution_sequence():
    """
    Based on qualitative analysis and simulation, this function returns the identified sequence.
    My detailed step-by-step reasoning pointed to a unique combination of (plot, variable, parameter).
    - Plots 3 & 4 show kinks characteristic of quarantine parameter changes (q_f, q_s).
    - Plot 6 shows an epidemic slowed down, consistent with a longer incubation period (a_i).
    - Plots 5, 7, 9 are increasing cumulative plots (Deaths, Costs, etc.).
        - Plot 5 (D vs f_s): Deaths are highly sensitive to the fraction developing severe symptoms.
        - Plot 7 (Cl vs c_l): The kink and scaling behavior match the cost of lost productivity.
        - Plot 9 (Ch vs c_h): A smooth S-curve representing cumulative hospital costs (a scaling effect).
    - Plots 1, 2, 8 are all Susceptible plots, distinguished by the impact of the varied parameter.
        - Plot 8 (S vs mu): The effect of baseline mortality is subtle, matching the closely-spaced curves.
        - Plot 2 (S vs mu_s): The mortality of severe cases has a surprisingly strong effect on the overall dynamics, leading to the dramatic drop.
        - Plot 1 (S vs beta_h): The effect of hospital transmission is significant, but less dramatic than that of mu_s.
    """
    
    # Mapping parameter names to their integer IDs
    param_map = {
        'mu': 1, 'mu_s': 2, 'mu_n': 3, 'a_i': 5, 'f_s': 6, 'c_l': 7, 
        'mu_h': 8, 'beta_h': 9, 'r_b': 11, 'c_h': 12, 'q_s': 13, 
        'q_l': 14, 'q_f': 15
    }

    # The identified parameters for plots 1 through 9
    plot_params = [
        'beta_h',  # Plot 1
        'mu_s',    # Plot 2
        'q_f',     # Plot 3
        'q_s',     # Plot 4
        'f_s',     # Plot 5
        'a_i',     # Plot 6
        'c_l',     # Plot 7
        'mu',      # Plot 8
        'c_h'      # Plot 9
    ]
    
    # Convert parameter names to their IDs
    solution_sequence = [param_map[p] for p in plot_params]
    
    return solution_sequence

# Get the final answer
final_sequence = get_solution_sequence()

# Print the result in the required format
# The final result is a sequence of numbers, representing the identified parameter ID for each plot.
print("Final Answer Sequence:")
# The format {p1, p2, ..., p9} is specified for the final output.
print(f"{{{', '.join(map(str, final_sequence))}}}")