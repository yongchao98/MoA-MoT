# Patient's blood pressure readings
systolic_rest = 135
diastolic_rest = 85
systolic_stand = 150
diastolic_stand = 90

# Calculate Pulse Pressure (PP)
pp_rest = systolic_rest - diastolic_rest
pp_stand = systolic_stand - diastolic_stand

# Calculate Mean Arterial Pressure (MAP)
map_rest = diastolic_rest + (pp_rest / 3)
map_stand = diastolic_stand + (pp_stand / 3)

print("Calculating the patient's Mean Arterial Pressure (MAP) at rest:")
print(f"MAP = Diastolic BP + (Systolic BP - Diastolic BP) / 3")
print(f"Resting MAP = {diastolic_rest} + ({systolic_rest} - {diastolic_rest}) / 3")
print(f"Resting MAP = {diastolic_rest} + {pp_rest} / 3")
print(f"Resting MAP = {map_rest:.2f} mmHg\n")

print("Calculating the patient's Mean Arterial Pressure (MAP) on standing:")
print(f"MAP = Diastolic BP + (Systolic BP - Diastolic BP) / 3")
print(f"Standing MAP = {diastolic_stand} + ({systolic_stand} - {diastolic_stand}) / 3")
print(f"Standing MAP = {diastolic_stand} + {pp_stand} / 3")
print(f"Standing MAP = {map_stand:.2f} mmHg")