import sys
# This script calculates the estimated age since eclosion for a mating male
# of Tridactylophagus tartari and Periplaneta americana.

# Step 1: Analyze Tridactylophagus tartari (a parasitoid wasp)
print("Analysis for Tridactylophagus tartari:")

# For a holometabolous insect (complete metamorphosis), "eclosion" is the emergence of the adult from the pupa.
# In this species, males emerge before females (protandry) and wait to mate.
male_emergence_lead_time_days = 3

print(f"This is a holometabolous insect, so eclosion is the emergence from the pupal case.")
print(f"Research shows that males emerge approximately {male_emergence_lead_time_days} days before females.")
print("They mate with females shortly after the females emerge. Therefore, the male's age since eclosion is his emergence lead time.")
print(f"Final Equation (T. tartari): Male Age = Lead Time = {male_emergence_lead_time_days} days")
age_t_tartari = f"{male_emergence_lead_time_days} days"
print("-" * 40)

# Step 2: Analyze Periplaneta americana (American cockroach)
print("\nAnalysis for Periplaneta americana:")

# For a hemimetabolous insect (incomplete metamorphosis), "eclosion" is hatching from the egg.
# The insect then undergoes a long nymphal stage before becoming a sexually mature adult.
avg_nymphal_period_months = 6

print("This is a hemimetabolous insect, so eclosion is hatching from the egg.")
print(f"The insect then develops as a nymph for an average of {avg_nymphal_period_months} months before molting into a sexually mature adult.")
print("Mating occurs after this final molt. Therefore, the male's age since eclosion is the length of its nymphal period.")
print(f"Final Equation (P. americana): Male Age = Nymphal Period = {avg_nymphal_period_months} months")
age_p_americana = f"{avg_nymphal_period_months} months"
print("-" * 40)

# Step 3: Combine the results to find the correct answer choice.
print("\nConclusion:")
print(f"The best average estimate of the male's age since eclosion is {age_t_tartari} for Tridactylophagus tartari and {age_p_americana} for Periplaneta americana.")