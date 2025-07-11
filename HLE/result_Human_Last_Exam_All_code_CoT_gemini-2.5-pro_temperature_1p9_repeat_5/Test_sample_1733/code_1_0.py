import math

# Patient Information
weight_kg = 12
dehydration_percentage = 0.10 # 10%

# --- 1. Calculate Initial Resuscitation Volume ---
bolus_rate_per_kg = 30
resuscitation_volume = weight_kg * bolus_rate_per_kg

print("1. Resuscitation Volume Calculation:")
print(f"   Volume = Body Weight (kg) * Bolus Rate (mL/kg)")
print(f"   Volume = {weight_kg} * {bolus_rate_per_kg} = {int(resuscitation_volume)} mL")
print("-" * 30)

# --- 2. Calculate Daily Maintenance Fluid Volume (Holliday-Segar) ---
maintenance_volume_24h = 0
remaining_weight = weight_kg
# First 10 kg
first_tier_weight = min(remaining_weight, 10)
first_tier_fluid = first_tier_weight * 100
maintenance_volume_24h += first_tier_fluid
remaining_weight -= first_tier_weight

# Next 10 kg (up to 20 kg total)
second_tier_weight = min(remaining_weight, 10)
second_tier_fluid = second_tier_weight * 50
maintenance_volume_24h += second_tier_fluid
remaining_weight -= second_tier_weight

# Above 20 kg
if remaining_weight > 0:
    third_tier_fluid = remaining_weight * 20
    maintenance_volume_24h += third_tier_fluid

print("2. Daily Maintenance Fluid Volume Calculation (Holliday-Segar):")
print(f"   For first {first_tier_weight} kg: {first_tier_weight} * 100 = {int(first_tier_fluid)} mL")
if second_tier_weight > 0:
    print(f"   For next {second_tier_weight} kg: {second_tier_weight} * 50 = {int(second_tier_fluid)} mL")
print(f"   Total Daily Maintenance = {int(first_tier_fluid)} + {int(second_tier_fluid)} = {int(maintenance_volume_24h)} mL")
print("-" * 30)


# --- 3. Calculate Total Deficit Replacement Fluid Volume ---
# 1 kg weight loss due to dehydration is equivalent to 1000 mL of fluid deficit.
deficit_volume = weight_kg * dehydration_percentage * 1000

print("3. Total Deficit Replacement Fluid Volume Calculation:")
print(f"   Volume = Body Weight (kg) * Dehydration % * 1000 mL/kg")
print(f"   Volume = {weight_kg} * {dehydration_percentage} * 1000 = {int(deficit_volume)} mL")
print("-" * 30)


# --- Final Answer ---
# Format: Resuscitation Volume, Daily Maintenance Volume, Total Deficit Volume
final_answer_str = f"{int(resuscitation_volume)},{int(maintenance_volume_24h)},{int(deficit_volume)}"
print("Final Answer (Resuscitation, Maintenance, Deficit):")
print(final_answer_str)
print(f"<<<{final_answer_str}>>>")