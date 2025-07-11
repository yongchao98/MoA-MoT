# The case describes a 43-year-old female with a 7-year history of symptoms
# whose condition resolved after 3 weeks before rapidly deteriorating.
patient_age = 43
history_duration_years = 7
symptom_resolution_weeks = 3

# This presentation is highly characteristic of a severe flare of Lupus Nephritis.
# To determine the immunological cause of such a flare, specific markers of
# SLE disease activity are required.
#
# Anti-dsDNA antibodies are specific to SLE, and their titers (levels) directly
# correlate with disease activity, particularly in cases of Lupus Nephritis.
# A sharp increase in this parameter would have been the strongest indicator
# of the cause behind the rapid renal decline.
most_indicative_lab_parameter = "Anti-double-stranded DNA (anti-dsDNA) antibodies"

# Print out the case numbers and the final answer.
print(f"For the patient who is {patient_age} years old with a {history_duration_years}-year history,")
print(f"and whose rebound symptoms resolved over {symptom_resolution_weeks} weeks, the following applies:")
print("\nThe lab parameter that could have best indicated the cause of the rapid renal function decline is:")
print(most_indicative_lab_parameter)
