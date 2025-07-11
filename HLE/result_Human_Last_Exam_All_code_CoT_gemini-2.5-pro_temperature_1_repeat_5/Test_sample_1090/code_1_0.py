# Patient's blood pressure values
resting_systolic_bp = 135
resting_diastolic_bp = 85
standing_systolic_bp = 150
standing_diastolic_bp = 90

# Calculate pulse pressure (Systolic BP - Diastolic BP)
resting_pulse_pressure = resting_systolic_bp - resting_diastolic_bp
standing_pulse_pressure = standing_systolic_bp - standing_diastolic_bp

# Print the results
print("Calculating the patient's pulse pressure:")
print(f"Resting Pulse Pressure = {resting_systolic_bp} mmHg - {resting_diastolic_bp} mmHg = {resting_pulse_pressure} mmHg")
print(f"Standing Pulse Pressure = {standing_systolic_bp} mmHg - {standing_diastolic_bp} mmHg = {standing_pulse_pressure} mmHg")
