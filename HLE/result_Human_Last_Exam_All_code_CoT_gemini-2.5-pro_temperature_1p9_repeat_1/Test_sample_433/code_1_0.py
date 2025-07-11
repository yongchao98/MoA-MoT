# Define the parameters based on the question
container_type = "single dose container ampule"
environment_condition = "sterile environment (ISO Class 5 or better)"
time_period_unit = "hours"

# According to USP General Chapter <797>, once a single-dose container
# is punctured or opened within an ISO Class 5 or better environment,
# it must be used within the following timeframe.
bud_in_hours = 12

# Print the explanation and the final answer
print(f"The Beyond-Use Date (BUD) for a {container_type} from the time of puncture in a {environment_condition} is {bud_in_hours} {time_period_unit}.")